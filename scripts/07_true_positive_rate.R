#-------------------------------------------------------------------------------
#         Process GWAS results from all methods: True positive rate
#-------------------------------------------------------------------------------
# we used an linkage disequilibrium (r^2) threshold of 0.1 
ld_threshold <- 0.1

# To avoid possible issues, it's good to close any GDS file opened before running the script
gdsfmt::showfile.gds(closeall = TRUE, verbose = F)

path = "./simulation/"
ls <- dir(path)
folder <- ls[!grepl("null",ls)]
lfolder <- length(folder)
results <- vector("list", lfolder)
# loop over all different folders containing results for the different settings
for (z in seq_along(folder)){
  species <- ifelse(grepl("maize",folder[z]), "maize", "soy")
  size <- ifelse(grepl("500",folder[z]), 500, 2815)
  # indicate the GDS file that will be opened for calculating LD
  if (species == "maize") {
    gds <- ifelse(grepl("500",folder[z]),
                  "./data/maize/ames_500_mac5_ld09_imputed.gds",
                  "./data/maize/ames_2815_mac5_ld09_imputed.gds")
  } else {
    gds <- ifelse(grepl("500",folder[z]),
                  "./data/soybean/soy_500_mac5_ld09.gds",
                  "./data/soybean/soy_2815_mac5_ld09.gds")
  }
  genofile <- snpgdsOpen(gds)
  # list all SNPs from the GDS file
  snps <- read.gdsn(index.gdsn(genofile, "snp.id"))
  scenario <- gsub("_.*", "",folder[z])
  # read the list of SNPs used as QTNs
  qtn <- fread(
    paste0(path, folder[z], "/Additive_QTNs.txt"),
    header = TRUE,
    sep = "\t",
    data.table = F
  )
  # list results from the multivariate stepwise model
  mtstep <- dir(paste0(path,folder[z], "/results"))
  # identify the result from the permutation
  permu <- dir(paste0(path,folder[z],"/results"))
  permu <- paste0(path,folder[z],"/results/",permu[grepl("_permutations.txt", permu)])
  pvalue <- unlist(fread(here(permu), data.table = F , nrows= 5))[-1:-4]
  mtstep <- mtstep[grepl("_Manova.txt",mtstep)]
  # list results from the univariate stepwise model
  ststep <- dir(paste0(path,folder[z], "/results_ST"))
  ststep <- ststep[!grepl("Marker|CI_|Permuted",ststep)]
  st <- gsub(".*trait|_Herit.*", "", ststep)
  st <- Reduce(rbind,strsplit(st,"Rep"))
  # list results from GEMMA
  gemma_results <- dir(paste0(path,folder[z], "/results_gemma/output"))
  gemma_results <- gemma_results[!grepl("log.txt",gemma_results)]
  
  n_gwas <- length(gemma_results)
  gwas_gemma <- vector("list", n_gwas)
  detection <- tibble(
    gemma_1 = logical(n_gwas),
    mt_1 = logical(n_gwas),
    t1 = logical(n_gwas),
    gemma_2 = logical(n_gwas),
    mt_2 = logical(n_gwas),
    t2 = logical(n_gwas),
    gemma_3 = logical(n_gwas),
    mt_3 = logical(n_gwas),
    t3 = logical(n_gwas))
  # loop over all the 100 settings
  for (n in 1:n_gwas){
    # get the replication number
    rep <- gsub("Rep|.assoc.txt|_Herit.*", "", gemma_results[n])
    # read the results from GEMMA
    gemma <- fread(
      paste0(path, folder[z], "/results_gemma/output/", gemma_results[n]),
      header = TRUE,
      sep = "\t",
      select = c("chr", "rs", "ps", "p_lrt"),
      data.table = F
    ) %>%
      arrange(p_lrt) %>%
      filter(p_lrt <= pvalue, !rs %in% qtn$snp)
    # read the results from the multivariate stepwise model
    mt <- fread(
      paste0(path, folder[z], "/results/", mtstep[n]),
      header = TRUE,
      sep = "\t",
      select = c("Chr", "SiteID", "Position", "probF"),
      data.table = F
    )
    # read univariate results for trait 1
    t1 <- fread(
      paste0(path, folder[z], "/results_st/", ststep[st[,2] %in% rep & st[,1] %in% 1]),
      header = TRUE,
      sep = "\t",
      data.table = F
    )
    # read univariate results for trait 2
    t2 <- fread(
      paste0(path, folder[z], "/results_st/", ststep[st[,2] %in% rep & st[,1] %in% 2]),
      header = TRUE,
      sep = "\t",
      data.table = F
    )
    # read univariate results for trait 3 when simulating partial pleiotropy
    if (scenario == "partially"){
      t3 <- fread(
        paste0(path, folder[z], "/results_st/", ststep[st[,2] %in% rep & st[,1] %in% 3]),
        header = TRUE,
        sep = "\t",
        data.table = F
      )
      t3 <- t3[grepl("Trait_3",t3$Trait) & !grepl("mean|Error",t3$Name),c("Locus", "Name", "Position", "pr>F" )]
      colnames(t3) = colnames(gemma)
    }
    # remove the mean and error component from the stepwise model result 
    t2 <- t2[grepl("Trait_2",t2$Trait) & !grepl("mean|Error",t2$Name),c("Locus", "Name", "Position", "pr>F" )]
    t1 <- t1[grepl("Trait_1",t1$Trait) & !grepl("mean|Error",t1$Name),c("Locus", "Name", "Position", "pr>F" )]
    colnames(t1) = colnames(t2) = colnames(mt) = colnames(gemma)
    # identifying true positives in the spurious pleiotropy scenario
    if (scenario == "LD") {
      detected_mt_1 <- c()
      detected_gemma_1 <- c()
      detected_t1 <- c()
      detected_mt_2 <- c()
      detected_gemma_2 <- c()
      detected_t2 <- c()
      l <- 1
      qtnlist <- split(qtn, qtn$trait)
      # loop over all the QTNs, if any of the SNPs detected as significant are 
      # above the LD threshold (r2 > 0.1) it will be called a true positive
      # loop over QTNs controlling trait 1
      for(i in qtnlist$trait_1$snp){
        #MT
        if (i %in% mt$rs) {
          detected_mt_1[l] <- TRUE
        } else {
          set <- snps[snps %in% c(i, mt$rs)] 
          # calculate all pairwise LD
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_mt_1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #GEMMA
        if (i %in% gemma$rs) {
          detected_gemma_1[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, gemma$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_gemma_1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T1
        if (i %in% t1$rs) {
          detected_t1[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t1$rs))] 
          
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      names(detected_mt_1) <- qtnlist$trait_1$snp
      names(detected_t1) <- qtnlist$trait_1$snp
      names(detected_gemma_1) <- qtnlist$trait_1$snp
      l <- 1
      # loop over QTNs controlling trait 2
      for(i in qtnlist$trait_2$snp){
        #MT
        if (i %in% mt$rs) {
          detected_mt_2[l] <- TRUE
        } else {
          set <- snps[snps %in% c(i, mt$rs)] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_mt_2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #GEMMA
        if (i %in% gemma$rs) {
          detected_gemma_2[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, gemma$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_gemma_2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T2
        if (i %in% t2$rs) {
          detected_t2[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t2$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      names(detected_mt_2) <- qtnlist$trait_2$snp
      names(detected_t2) <- qtnlist$trait_2$snp
      names(detected_gemma_2) <- qtnlist$trait_2$snp
      # create a data frame when n = 1 and save it as "res". All other replicates
      # are added as additional columns (temp)
      if(n==1) {
        res <- rbind(
          data.frame(trait = 1, snp = c(1, 2, 3), maf = qtnlist$trait_1$maf, AF = qtnlist$trait_1$additive_effect,detected_mt_1, detected_t1, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 2, snp = c(1, 2, 3), maf = qtnlist$trait_2$maf, AF = qtnlist$trait_2$additive_effect,detected_mt_2, detected_t2, detected_gemma_2) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n))
        )
      } else {
        temp <- rbind(
          data.frame(trait = 1, snp = c(1, 2, 3), maf = qtnlist$trait_1$maf, AF = qtnlist$trait_1$additive_effect,detected_mt_1, detected_t1, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 2, snp = c(1, 2, 3), maf = qtnlist$trait_2$maf, AF = qtnlist$trait_2$additive_effect,detected_mt_2, detected_t2, detected_gemma_2) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n))
        )
        res <- res %>% left_join(temp,by = c("trait", "snp", "maf", "AF", "model") )
      }
      # identifying true positives in the partial pleiotropy scenario
    } else if (scenario == "partially") {
      detected_mt_1 <- c()
      detected_gemma_1 <- c()
      detected_t1 <- c()
      detected_mt_2 <- c()
      detected_gemma_2 <- c()
      detected_t2 <- c()
      detected_mt_3 <- c()
      detected_gemma_3 <- c()
      detected_t3 <- c()
      l <- 1
      qtnlist <- split(qtn, qtn$trait)
      for(i in qtnlist$trait_1$snp){
        #MT
        if (i %in% mt$rs) {
          detected_mt_1[l] <- TRUE
        } else {
          set <- snps[snps %in% c(i, mt$rs)] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_mt_1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #GEMMA
        if (i %in% gemma$rs) {
          detected_gemma_1[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, gemma$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_gemma_1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T1
        if (i %in% t1$rs) {
          detected_t1[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t1$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD  ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      names(detected_mt_1) <- qtnlist$trait_1$snp
      names(detected_t1) <- qtnlist$trait_1$snp
      names(detected_gemma_1) <- qtnlist$trait_1$snp
      l <- 1
      for(i in qtnlist$trait_2$snp){
        #MT
        if (i %in% mt$rs) {
          detected_mt_2[l] <- TRUE
        } else {
          set <- snps[snps %in% c(i, mt$rs)] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_mt_2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #GEMMA
        if (i %in% gemma$rs) {
          detected_gemma_2[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, gemma$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_gemma_2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T2
        if (i %in% t2$rs) {
          detected_t2[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t2$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      names(detected_mt_2) <- qtnlist$trait_2$snp
      names(detected_t2) <- qtnlist$trait_2$snp
      names(detected_gemma_2) <- qtnlist$trait_2$snp
      l <- 1
      for(i in qtnlist$trait_3$snp){
        #MT
        if (i %in% mt$rs) {
          detected_mt_3[l] <- TRUE
        } else {
          set <- snps[snps %in% c(i, mt$rs)] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_mt_3[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #GEMMA
        if (i %in% gemma$rs) {
          detected_gemma_3[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, gemma$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_gemma_3[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T3
        if (i %in% t3$rs) {
          detected_t3[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t3$rs))]
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t3[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      names(detected_mt_3) <- qtnlist$trait_2$snp
      names(detected_t3) <- qtnlist$trait_2$snp
      names(detected_gemma_3) <- qtnlist$trait_2$snp
      if(n==1) {
        res <- rbind(
          data.frame(trait = 1, snp = c(1, 2, 3), maf = qtnlist$trait_1$maf, AF = qtnlist$trait_1$additive_effect,detected_mt_1, detected_t1, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 2, snp = c(1, 2, 3), maf = qtnlist$trait_2$maf, AF = qtnlist$trait_2$additive_effect,detected_mt_2, detected_t2, detected_gemma_2) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 3, snp = c(1, 2, 3), maf = qtnlist$trait_3$maf, AF = qtnlist$trait_3$additive_effect,detected_mt_3, detected_t3, detected_gemma_3) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n))
        )
      } else {
        temp <- rbind(
          data.frame(trait = 1, snp = c(1, 2, 3), maf = qtnlist$trait_1$maf, AF = qtnlist$trait_1$additive_effect,detected_mt_1, detected_t1, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 2, snp = c(1, 2, 3), maf = qtnlist$trait_2$maf, AF = qtnlist$trait_2$additive_effect,detected_mt_2, detected_t2, detected_gemma_2) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 3, snp = c(1, 2, 3), maf = qtnlist$trait_3$maf, AF = qtnlist$trait_3$additive_effect,detected_mt_3, detected_t3, detected_gemma_3) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n))
        )
        res <- res %>% left_join(temp,by = c("trait", "snp", "maf", "AF", "model") )
      }
      # identifying true positives in the pleiotropy scenario
    } else {
      detected_mt_1 <- c()
      detected_gemma_1 <- c()
      detected_t1 <- c()
      detected_t2 <- c()
      l <- 1
      for(i in qtn$snp){
        #MT
        if (i %in% mt$rs) {
          detected_mt_1[l] <- TRUE
        } else {
          set <- snps[snps %in% c(i, mt$rs)] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_mt_1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #GEMMA
        if (i %in% gemma$rs) {
          detected_gemma_1[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, gemma$rs))] 
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_gemma_1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T1
        if (i %in% t1$rs) {
          detected_t1[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t1$rs))]
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t1[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        #T2
        if (i %in% t2$rs) {
          detected_t2[l] <- TRUE
        } else {
          set <- snps[which(snps %in% c(i, t2$rs))]
          ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
          colnames(ld) = rownames(ld) = set
          detected_t2[l] <- any(ld[i,-which(colnames(ld) == i)] > ld_threshold)
        }
        l <- l + 1
      }
      names(detected_mt_1) <- qtn$snp
      names(detected_t1) <- qtn$snp
      names(detected_t2) <- qtn$snp
      names(detected_gemma_1) <- qtn$snp
      
      if(n==1) {
        res <- rbind(
          data.frame(trait = 1, snp = c(1, 2, 3), maf = qtn$maf, AF = qtn$add_eff_t1,detected_mt_1, detected_t1, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 2, snp = c(1, 2, 3), maf = qtn$maf, AF = qtn$add_eff_t2,detected_mt_1, detected_t2, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n))
        )
      } else {
        temp <- rbind(
          data.frame(trait = 1, snp = c(1, 2, 3), maf = qtn$maf, AF = qtn$add_eff_t1,detected_mt_1, detected_t1, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n)),
          data.frame(trait = 2, snp = c(1, 2, 3), maf = qtn$maf, AF = qtn$add_eff_t2,detected_mt_1, detected_t2, detected_gemma_1) %>%
            pivot_longer(!snp & !trait & !AF & !maf, names_to = "model", values_to = paste0("count", n))
        )
        res <- res %>% left_join(temp,by = c("trait", "snp", "maf", "AF", "model") )
      }
    }
    print(n)
  }
  # summarize all results and save in a list
  results[[z]] <- res %>% mutate(sum = rowSums(across(where(is.logical))), species = species, size = size, architecture = scenario) %>% select(architecture, species, size, trait, snp, model,maf, AF, sum) 
  # close the GDS file
  gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
  cat("-----------done with ---------", z)
}
results <- Reduce(rbind, results)
# save the true positive rate
fwrite(results, "./results/TP.txt", sep = "\t", quote = F, row.names = F)
