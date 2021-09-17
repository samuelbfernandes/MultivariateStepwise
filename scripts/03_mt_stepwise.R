#-------------------------------------------------------------------------------
#         calling the multivariate multi-locus stepwise model from TASSEL
#-------------------------------------------------------------------------------
#--------------- multivariate analysis non-NULL settings -----------------------
#list all directories where phenotypes are saved
dir_path <- paste0(gsub(".*simulation", "./simulation", dir_ls(here("simulation"))), "/Phenotypes")
dir_path <- dir_path[!grepl("null", dir_path)]
#Loop over all non-null settings
for (d in 1:length(dir_path)) {
  #list all files (replicate) in a directory
  files <- dir_ls(dir_path[d])
  #the first file is used for permutation
  out <- paste0(dirname(dir_path[d]),"/results/", gsub(".*__|_H.*","",files[1]))
  dir_create(dirname(out))
  #get the plink bed files with QTNs removed
  bed <- dir_ls(dirname(dir_path[d]))
  bed <- gsub(".bed", "", bed[grepl("bed", bed)])
  #converting bed files to vcf to be used in TASSEL
  system(command = paste(path_to_plink, " --bfile", bed, "--recode vcf --out", bed), intern = F, ignore.stdout = T, ignore.stderr = T)
  vcf <- paste0(bed, ".vcf")
  #calling the MultivariateStepwisePlugin from TASSEL with 100 permutations
  system(command =
           paste(path_to_tassel,
                 "-fork1 -vcf", vcf,
                 "-fork2 -r ", files[1],
                 " -combine3 -input1 -input2 -intersect",
                 "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05",
                 "-saveToFile true -savePath", out,
                 "-parallel true -threads", ncores, "-endPlugin -export -runfork1 -runfork2"), intern = F)
  #extract the 5th pvalue (for alpha 0.05) from permutation pvalues
  pvalue <- unlist(fread(paste0(out, "_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
  #loop over the next 99 replicates using the above pvalue as entreLimit
  for (i in 2:100){
    out <- paste0(dirname(dir_path[d]),"/results/", gsub(".*__|_H.*","",files[i]))
    system(command = 
             paste(path_to_tassel, 
                   "-fork1 -vcf", vcf,
                   "-fork2 -r ", files[i],
                   " -combine3 -input1 -input2 -intersect", 
                   "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                   "-saveToFile true -savePath", out,
                   "-parallel true -threads",  ncores, " -endPlugin",
                   "-export -runfork1 -runfork2"), intern = F)
  }
}

#--------------- multivariate analysis NULL settings ----------------------------
#---- creating folders for storing all NULL setting results ----
dir_create(here("simulation/null_cor_M/results"))
dir_create(here("simulation/null_cor_S/results"))
dir_create(here("simulation/null_M/results"))
dir_create(here("simulation/null_S/results"))
dir_create(here("simulation/null_cor_M_2815/results"))
dir_create(here("simulation/null_cor_S_2815/results"))
dir_create(here("simulation/null_M_2815/results"))
dir_create(here("simulation/null_S_2815/results"))

#------ Maize correlated phenotypes sample size 500 ------
files_cm <- dir_ls(here("simulation/null_cor_M/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/maize/ames_500_mac5_ld09_imputed.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_cm[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", paste0("./simulation/null_cor_M/results/", gsub(".*__|_H.*","",files_cm[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_cor_M/results/Rep100_permutations.txt"), data.table = F , nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/maize/ames_500_mac5_ld09_imputed.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_cm[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", paste0("./simulation/null_cor_M/results/", gsub(".*__|_H.*","",files_cm[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Soybean correlated phenotypes sample size 500 ------
files_cs <- dir_ls(here("simulation/null_cor_S/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/soybean/soy_500_mac5_ld09.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_cs[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", 
               paste0("./simulation/null_cor_S/results/", gsub(".*__|_H.*","",files_cs[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_cor_S/results/Rep100_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/soybean/soy_500_mac5_ld09.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_cs[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", 
                 paste0("./simulation/null_cor_S/results/", gsub(".*__|_H.*","",files_cs[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Maize Uncorrelated phenotypes sample size 500 ------
files_m <- dir_ls(here("simulation/null_M/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/maize/ames_500_mac5_ld09_imputed.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_m[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", paste0("./simulation/null_M/results/",
                                                    gsub(".*__|_H.*","",files_m[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_M/results/Rep100_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/maize/ames_500_mac5_ld09_imputed.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_m[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", paste0("./simulation/null_M/results/", gsub(".*__|_H.*","",files_m[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Soybean Uncorrelated phenotypes sample size 500 ------
files_s <- dir_ls(here("simulation/null_S/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/soybean/soy_500_mac5_ld09.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_s[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", 
               paste0("./simulation/null_S/results/", gsub(".*__|_H.*","",files_s[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_S/results/Rep100_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/soybean/soy_500_mac5_ld09.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_s[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", 
                 paste0("./simulation/null_S/results/", gsub(".*__|_H.*","",files_s[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Maize correlated phenotypes sample size 2815 ------
files_cm <- dir_ls(here("simulation/null_cor_M_2815/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/maize/ames_2815_mac5_ld09_imputed.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_cm[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", paste0("./simulation/null_cor_M_2815/results/",
                                                    gsub(".*__|_H.*","",files_cm[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_cor_M_2815/results/Rep100_permutations.txt"), data.table = F , nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/maize/ames_2815_mac5_ld09_imputed.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_cm[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", paste0("./simulation/null_cor_M_2815/results/",
                                                      gsub(".*__|_H.*","",files_cm[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Soybean correlated phenotypes  sample size 2815 ------
files_cs <- dir_ls(here("simulation/null_cor_S_2815/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/soybean/soy_2815_mac5_ld09.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_cs[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", 
               paste0("./simulation/null_cor_S_2815/results/", gsub(".*__|_H.*","",files_cs[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_cor_S_2815/results/Rep100_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/soybean/soy_2815_mac5_ld09.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_cs[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", 
                 paste0("./simulation/null_cor_S_2815/results/", gsub(".*__|_H.*","",files_cs[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Maize Uncorrelated phenotypes  sample size 2815 ------
files_m <- dir_ls(here("simulation/null_M_2815/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/maize/ames_2815_mac5_ld09_imputed.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_m[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", paste0("./simulation/null_M_2815/results/", gsub(".*__|_H.*","",files_m[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_M_2815/results/Rep100_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/maize/ames_2815_mac5_ld09_imputed.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_m[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", paste0("./simulation/null_M_2815/results/", gsub(".*__|_H.*","",files_m[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

#------ Soybean Uncorrelated phenotypes  sample size 2815 ------
files_s <- dir_ls(here("simulation/null_S_2815/Phenotypes"))
system(command = 
         paste(path_to_tassel,
               "-fork1 -h ./data/soybean/soy_2815_mac5_ld09.hmp.txt",
               "-fork2 -r ", gsub(".*simulation", "./simulation", files_s[1]),
               " -combine3 -input1 -input2 -intersect", 
               "-MultivariateStepwisePlugin -usePerm true -nPerm 100 -permAlpha 0.05", 
               "-saveToFile true -savePath", 
               paste0("./simulation/null_S_2815/results/", gsub(".*__|_H.*","",files_s[1])),
               "-parallel true -threads ", ncores, " -endPlugin",
               "-export -runfork1 -runfork2"), intern = F)

pvalue <- unlist(fread(here("simulation/null_S_2815/results/Rep100_permutations.txt"), data.table = F, nrows= 5))[-1:-4]
for (i in 2:100){
  system(command = 
           paste(path_to_tassel,
                 "-fork1 -h ./data/soybean/soy_2815_mac5_ld09.hmp.txt", 
                 "-fork2 -r ", gsub(".*simulation", "./simulation", files_s[i]),
                 " -combine3 -input1 -input2 -intersect", 
                 "-MultivariateStepwisePlugin -usePerm false -enterLimit", pvalue, "-exitLimit", pvalue * 2, 
                 "-saveToFile true -savePath", 
                 paste0("./simulation/null_S_2815/results/", gsub(".*__|_H.*","",files_s[i])),
                 "-parallel true -threads ", ncores, " -endPlugin",
                 "-export -runfork1 -runfork2"), intern = F)
}

