#-------------------------------------------------------------------------------
#         Auxiliar functions to run GEMMA and to obtain a kinship Matrix
#-------------------------------------------------------------------------------

# Function to call GEMMA from R
gemma <- function(gemma_path = NULL,
                  path_out = getwd(),
                  geno = NULL,
                  fixed = NULL,
                  kin = NULL,
                  d = NULL,
                  u = NULL,
                  lmm = 2,
                  miss = 0.05,
                  maf = 0.01,
                  r2 = 0.9999,
                  trait_col = NULL,
                  out_name = NULL,
                  verbose = T){
  home_path <- getwd()
  setwd(path_out)
  if (!is.null(kin)) {
    d <- NULL
    u <- NULL}
  if (is.null(fixed)) {
    if (!is.null(d) &
        !is.null(u)) {
      system(command = paste(gemma_path, 
                             "--bfile",
                             geno,
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-n", paste(trait_col, collapse = " "),
                             "-d", d,
                             "-u", u,
                             "-o", out_name),
             ignore.stdout = !verbose,
             ignore.stderr = !verbose)
    } else {
      system(command = paste(gemma_path, 
                             "--bfile",
                             geno,
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-n", paste(trait_col, collapse = " "),
                             "-k", kin,
                             "-o", out_name),
             ignore.stdout = !verbose,
             ignore.stderr = !verbose)
    }
  } else {
    if (!is.null(d) &
        !is.null(u)) {
      system(command = paste(gemma_path, 
                             "--bfile",
                             geno,
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-c",fixed,
                             "-n", paste(trait_col, collapse = " "),
                             "-d", d,
                             "-u", u,
                             "-o", out_name),
             ignore.stdout = !verbose,
             ignore.stderr = !verbose)
    } else {
      system(command = paste(gemma_path, 
                             "--bfile",
                             geno,
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-c",fixed,
                             "-n", paste(trait_col, collapse = " "),
                             "-k", kin,
                             "-o", out_name),
             ignore.stdout = !verbose,
             ignore.stderr = !verbose)
    }
  }
  setwd(home_path)
}

# Function to obtain a kinship matrix from GEMMA
relatedness <- function(path_to_tassel,
                        path_to_plink,
                        path_to_gemma,
                        data){
  # convert from HapMap to Plink Bed files
  system(command = paste(path_to_tassel, "-fork1 -h", paste0(data,".hmp.txt"), "-export ", paste0(data,".vcf"), "-exportType VCF -runfork1"), intern = F)
  system(command = paste(path_to_plink, "--vcf", paste0(data,".vcf"), "--make-bed --out", data), intern = F, ignore.stdout = T, ignore.stderr = T)
  df <- fread(paste0(data,".fam"), data.table = F)
  df[,6] <- 1
  fwrite(df, paste0(data,".fam"), quote = F, row.names = F, sep = "\t", col.names = F)
  setwd(dirname(data))
  file <- gsub(".*/", "", data)
  # relatedness and eigen decomposition
  system(command = paste(path_to_gemma, "-bfile", file, "-gk 1 -r2 1 -maf 0 -miss 1 -o", file), intern = F, ignore.stdout = T, ignore.stderr = T)
  system(command = paste(path_to_gemma, "-bfile", file, "-k", paste0("./output/",file,".cXX.txt"),
                         "-eigen -r2 0.9999999 -maf 0.0001 -miss 0.9999 -o", file), intern = F, ignore.stdout = T, ignore.stderr = T)
  file_move(paste0("./output/", file, ".eigenD.txt"), paste0("./", file, ".eigenD.txt"))
  file_move(paste0("./output/", file, ".eigenU.txt"), paste0("./", file, ".eigenU.txt"))
  unlink("output", recursive = T)
  setwd(here())
}
