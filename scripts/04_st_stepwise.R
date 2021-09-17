#-------------------------------------------------------------------------------
#         calling the Univariate multi-locus stepwise model from TASSEL
#-------------------------------------------------------------------------------
#list all directories where phenotypes are saved
dir_path <- paste0(gsub(".*simulation", "./simulation", dir_ls(here("simulation"))), "/Phenotypes")
dir_path <- dir_path[!grepl("null", dir_path)]
#Loop over all non-null settings
for (d in 1:length(dir_path)) {
  #list all files (replicate) in a directory
  files <- dir_ls(dir_path[d])
  dir_create(paste0(dirname(dir_path[d]), "/results_ST"))
  size <- gsub(".*_|/.*","",dir_path[d])
  scenario <- gsub("_.*|./simulation/","",dir_path[d])
  species <- gsub(".*_","",gsub(paste0("_",size, ".*"),"",dir_path[d]))
  vcf <- dir_ls(dirname(dir_path[d]))
  vcf <- vcf[grepl("vcf", vcf)]
  
  # permutation p-values for the univariate model were obtained on TASSEL GUI
  # the original implementation only outputs the permutation for the last trait
  # to get the p-value for all traits we need to look at the enter limit printed
  # on the log file.
  pvalue <- fread(here("./ST_permutation_pvalues.txt"), data.table = F)
  pvalue <- pvalue[pvalue$setting %in% gsub("./simulation/|/Phenotypes","",dir_path[d]),]
  
  # make the output directory the working directory 
  setwd(paste0(dirname(dir_path[d]), "/results_ST"))
  # make a list of files to parallelize using mclapply
  files_fam <- split(files, 1:length(files))
  # run the univariate stepwise model on each simulated experiment in parallel
  mclapply(files_fam, function(x) {
    pheno <- fread(here(x), data.table = F)
    # the scenario of partial pleiotropy contains three traits. Thus, it calls 
    # TASSEL 3x, one for each trait.
    if (scenario == "partially") {
      t1 <- paste0("trait1", gsub(".*__", "",x))
      t2 <- paste0("trait2", gsub(".*__", "",x))
      t3 <- paste0("trait3", gsub(".*__", "",x))
      # simplePHENOTYPES outputs all traits from a given experiment in a single file
      # to run the univariate stepwise model with a different permutation threshold
      # each trait needs to be saved in a separate file
      fwrite(pheno[, 1:2], t1, sep = "\t", quote = F, row.names = F)
      fwrite(pheno[, c(1,3)], t2, sep = "\t", quote = F, row.names = F)
      fwrite(pheno[, c(1,4)], t3, sep = "\t", quote = F, row.names = F)
      system(command = 
               paste(path_to_tassel, 
                     "-fork1 -vcf", paste0("../../../",vcf),
                     "-fork2 -r ", t1,
                     " -combine3 -input1 -input2 -intersect", 
                     "-StepwiseOLSModelFitterPlugin -enter", pvalue$trait1 , "-exit", pvalue$trait1 *2, 
                     "-endPlugin",
                     "-export -runfork1 -runfork2"), intern = F) 
      system(command = 
               paste(path_to_tassel,
                     "-fork1 -vcf", paste0("../../../",vcf),
                     "-fork2 -r", t2,
                     " -combine3 -input1 -input2 -intersect", 
                     "-StepwiseOLSModelFitterPlugin -enter", pvalue$trait2 , "-exit", pvalue$trait2 *2, 
                     "-endPlugin",
                     "-export -runfork1 -runfork2"), intern = F) 
      system(command = 
               paste(path_to_tassel,
                     "-fork1 -vcf", paste0("../../../",vcf),
                     "-fork2 -r", t3,
                     " -combine3 -input1 -input2 -intersect", 
                     "-StepwiseOLSModelFitterPlugin -enter", pvalue$trait3 , "-exit", pvalue$trait3 *2, 
                     "-endPlugin",
                     "-export -runfork1 -runfork2"), intern = F) 
      unlink(c(t1,t2,t3), recursive = T, force = T)
    } else {
      t1 <- paste0("trait1", gsub(".*__", "",x))
      t2 <- paste0("trait2", gsub(".*__", "",x))
      fwrite(pheno[, 1:2], t1, sep = "\t", quote = F, row.names = F)
      fwrite(pheno[, c(1,3)], t2, sep = "\t", quote = F, row.names = F)
      system(command = 
               paste(path_to_tassel, 
                     "-fork1 -vcf", paste0("../../../",vcf),
                     "-fork2 -r ", t1,
                     " -combine3 -input1 -input2 -intersect", 
                     "-StepwiseOLSModelFitterPlugin -enter", pvalue$trait1 , "-exit", pvalue$trait1 *2, 
                     "-endPlugin",
                     "-export -runfork1 -runfork2"), intern = F) 
      system(command = 
               paste(path_to_tassel, 
                     "-fork1 -vcf", paste0("../../../",vcf),
                     "-fork2 -r", t2,
                     " -combine3 -input1 -input2 -intersect", 
                     "-StepwiseOLSModelFitterPlugin -enter", pvalue$trait2 , "-exit", pvalue$trait2 *2, 
                     "-endPlugin",
                     "-export -runfork1 -runfork2"), intern = F) 
      unlink(c(t1,t2), recursive = T, force = T)
      }
  },
  mc.preschedule = FALSE,
  mc.cores = ncores)
  setwd(here())
}

