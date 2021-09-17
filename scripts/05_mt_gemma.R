#-------------------------------------------------------------------------------
#         calling the Multivariate single locus model (mvLMM) from GEMMA
#-------------------------------------------------------------------------------
#list all directories where phenotypes are saved
dir_path <- paste0(gsub(".*simulation", "./simulation", dir_ls(here("simulation"))), "/Phenotypes")
dir_path <- dir_path[!grepl("null", dir_path)]
#Loop over all non-null settings
for (d in 1:length(dir_path)) {
  files <- dir_ls(dir_path[d])
  out <- paste0(dirname(dir_path[d]),"/results/", gsub(".*__|_H.*","",files[1]))
  dir_create(paste0(dirname(dir_path[d]), "/results_gemma"))
  size <- gsub(".*_|/.*","",dir_path[d])
  species <- gsub(".*_","",gsub(paste0("_",size, ".*"),"",dir_path[d]))
  # make the output directory the working directory 
  setwd(paste0(dirname(dir_path[d]), "/results_gemma"))
  # copying the eigen values and engen vectors of the kinship matrix of
  # maize and soybean to the respective working directory
  if (species == "maize") {
    bed <- paste0("../ames_", size, "_mac5_ld09_imputed_noQTN")
    file_copy(here(paste0("./data/maize/ames_", size, "_mac5_ld09_imputed.eigenD.txt")), "kinD.txt")
    file_copy(here(paste0("./data/maize/ames_", size, "_mac5_ld09_imputed.eigenU.txt")), "kinU.txt")
  } else {
    bed <- paste0("../soy_", size, "_mac5_ld09_noQTN")
    file_copy(here(paste0("./data/soybean/soy_", size, "_mac5_ld09.eigenD.txt")), "kinD.txt")
    file_copy(here(paste0("./data/soybean/soy_", size, "_mac5_ld09.eigenU.txt")), "kinU.txt")
  }
  # make a list of files to parallelize using mclapply
  files_fam <- split(files, 1:length(files))
  # run the univariate stepwise model on each simulated experiment in parallel
  mclapply(files_fam, function(x) {
    # converting and saving the phenotypic file to a ".fam" file to match with
    # the genomic data used by GEMMA
    pheno <- fread(here(x), data.table = F)
    fam <- fread(paste0(bed, ".fam"), data.table = F)
    fam <- merge(fam[,-1], pheno, by.x = "V2", by.y = "<Trait>", sort = F)
    fam <- fam[, c(1, 1, 2, 3, 4, 6:ncol(fam))]
    fwrite(fam, gsub(dirname(x), ".", gsub(".txt", ".fam" ,x)), sep = "\t", quote = T, col.names = F, row.names = F)
    # renaming BED files for each replicate
    file_copy(paste0(bed, ".bed"), gsub(dirname(x), ".", gsub(".txt", ".bed" ,x)))
    file_copy(paste0(bed, ".bim"), gsub(dirname(x), ".", gsub(".txt", ".bim" ,x)))
    # setting the number of traits to run in the multivariate GWAS model
    trait_col <- if (grepl("partially", out)) {
      c(1, 2, 3)
    } else{
      c(1, 2)
    }
    # calling GEMMA from R (see function gemma() inside "01_gemma.R" )
    gemma(path_to_gemma, 
          geno = gsub(".*Phenotypes/|.txt","",x),
          d = "kinD.txt",
          u = "kinU.txt",
          trait_col = trait_col,
          miss = 0.001,
          maf = 0.001,
          r2 = 0.999999,
          out_name = gsub(".*__|_Herit_0.3_0.8.txt", "", x)
    )
    unlink(c(
      paste0(gsub(".*Phenotypes/|.txt","",x), ".bed"),
      paste0(gsub(".*Phenotypes/|.txt","",x), ".bim"),
      paste0(gsub(".*Phenotypes/|.txt","",x), ".fam")
    ),
    recursive = T,
    force = T)
    print("Done! Starting next file...")
  },
  mc.preschedule = FALSE,
  mc.cores = ncores)
  setwd(here())
}
