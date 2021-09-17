#-------------------------------------------------------------------------------
#         Process GWAS results from all methods: False positive rate
#-------------------------------------------------------------------------------

#------------ Processing the results from all the NULL settings ----------------
path = "./simulation/"
ls <- dir(path)

folder <- ls[grepl("null",ls)]
lfolder <- length(folder)
results <- vector("list", lfolder)
# loop over all different folders containing results for the different settings
for (z in seq_along(folder)){
  species <- ifelse(grepl("M",folder[z]), "maize", "soy")
  size <- ifelse(grepl("2815",folder[z]), 2815, 500)
  scenario <- gsub("_.*", "",folder[z])
  # results from the multivariate stepwise model
  mtstep <- dir(paste0(path,folder[z], "/results"))
  mtstep <- mtstep[grepl("_Manova.txt",mtstep)]
  # the loop below will extract results from all 100 replicates simulated
  n_gwas <- 100
  FP_null <- tibble(
    species = species,
    size = size,
    cor = grepl("cor",folder[z]),
    mt = logical(n_gwas))
  for (n in 1:n_gwas){
    # read the results obtained with the multivariate stepwise model for 
    # the setting z and the replicate n
    mt <- fread(
      paste0(path, folder[z], "/results/", mtstep[n]),
      header = TRUE,
      sep = "\t",
      select = c("Chr", "SiteID", "Position", "probF"),
      data.table = F
    )
    colnames(mt) = c("Chr", "SiteID", "Position", "probF")
    FP_null$mt[n] <- (nrow(mt) > 0)
    print(n)
  }
  results[[z]] <- FP_null %>% group_by(species, size, cor) %>% summarise(sum = sum(mt)) %>% ungroup()
  cat("-----------done with ---------", z)
}
# convert the list into a data.frame
results <- Reduce(rbind, results)
# save the results from the null setting into a file
fwrite(results, "./results/FPR_NULL.txt", sep = "\t", quote = F, row.names = F)

#------------ Processing the results from all the non-NULL settings ----------------
# we used an linkage disequilibrium (r^2) threshold of 0.1 
ld_threshold <- 0.1

# To avoid possible issues, it's good to close any GDS file opened before running the script
gdsfmt::showfile.gds(closeall = TRUE, verbose = F)

folder <- ls[!grepl("null",ls)]
lfolder <- length(folder)
results <- vector("list", lfolder)
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
  mtstep <- dir(paste0(path,folder[z], "/results"))
  # identify the result from the permutation
  permu <- dir(paste0(path,folder[z], "/results"))
  permu <- paste0(path,folder[z], "/results/", permu[grepl("_permutations.txt", permu)])
  # select the 5th permutation pvalue to use as threshold for GEMMA
  pvalue <- unlist(fread(here(permu), data.table = F , nrows= 5))[-1:-4]
  mtstep <- mtstep[grepl("_Manova.txt", mtstep)]
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
  FP <- tibble(
    gemma = logical(n_gwas),
    mt = logical(n_gwas),
    t1 = logical(n_gwas),
    t2 = logical(n_gwas),
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
    # identifying false positives in the spurious pleiotropy scenario
    if (scenario == "LD") {
      FP_mt <- c()
      FP_gemma <- c()
      FP_t1 <- c()
      FP_t2 <- c()
      l <- 1
      # loop over all the SNPs detected as significant in a given replicate
      # call it false positive if it is not in LD (r2 > 0.1) with any of the QTNs
      for(i in mt$rs){
        #MT
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        # calculate all pairwise LD
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        # if the LD between a given SNP and all the QTNs is < 0.1 it will be called a false positive
        FP_mt[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in gemma$rs){
        #GEMMA
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_gemma[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t1$rs){
        #T1
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_t1[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t2$rs){
        #T2
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2 
        colnames(ld) = rownames(ld) = set
        FP_t2[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      # if any false positive was detected in a replicate, it is counted only once
      FP$gemma[n] <- any(FP_gemma)
      FP$mt[n] <- any(FP_mt)
      FP$t1[n] <- any(FP_t1)
      FP$t2[n] <- any(FP_t2)
      
      # identifying false positives in the partial pleiotropy scenario
    } else if (scenario == "partially") {
      FP_mt <- c()
      FP_gemma <- c()
      FP_t1 <- c()
      FP_t2 <- c()
      FP_t3 <- c()
      l <- 1
      for(i in mt$rs){
        #MT
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_mt[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in gemma$rs){
        #GEMMA
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_gemma[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t1$rs){
        #T1
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_t1[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t2$rs){
        #T2
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_t2[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t3$rs){
        #T3
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_t3[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      FP$gemma[n] <- any(FP_gemma)
      FP$mt[n] <- any(FP_mt)
      FP$t1[n] <- any(FP_t1)
      FP$t2[n] <- any(FP_t2)
      FP$t3[n] <- any(FP_t3)
      
      # identifying false positives in the pleiotropy scenario
    } else {
      FP_mt <- c()
      FP_gemma <- c()
      FP_t1 <- c()
      FP_t2 <- c()
      l <- 1
      for(i in mt$rs){
        #MT
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD  ^ 2
        colnames(ld) = rownames(ld) = set
        FP_mt[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in gemma$rs){
        #GEMMA
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_gemma[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t1$rs){
        #T1
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_t1[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      l <- 1
      for(i in t2$rs){
        #T2
        set <- snps[which(snps %in% c(i, qtn$snp))] 
        ld <- snpgdsLDMat(genofile, snp.id = set,  slide = -1, method = "r", verbose=F)$LD ^ 2
        colnames(ld) = rownames(ld) = set
        FP_t2[l] <- all(ld[i,-which(colnames(ld) == i)] < ld_threshold)
        l = l+1
      }
      FP$gemma[n] <- any(FP_gemma)
      FP$mt[n] <- any(FP_mt)
      FP$t1[n] <- any(FP_t1)
      FP$t2[n] <- any(FP_t2)
    }
    print(n)
  }
  # count the number of false positives on 100 replicates
  total_FP <- t(apply(FP, 2, sum))
  results[[z]] <- data.frame(species = species, size = size, architecture = scenario,  total_FP)
  # close the GDS file
  gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
  cat("-----------done with ---------", z)
}
results <- Reduce(rbind, results)
# save the non-null false positive rate
fwrite(results, "./results/FPR.txt", sep = "\t", quote = F, row.names = F)
