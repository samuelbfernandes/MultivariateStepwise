#-------------------------------------------------------------------------------
#         Example of a two step approach with kinship as a random effect
#-------------------------------------------------------------------------------
pheno <- "./simulation/LD_maize_500/Phenotypes/Simulated_Data__Rep100_Herit_0.3_0.8.txt"
hmp <- "./data/maize/ames_500_mac5_ld09_imputed.hmp.txt"

dir_create(here("two_step"))

# obtaining a kinship matrix. We are using TASSEL for convenience.
# However, if saved in TASSEL's square matrix format, any other kinship matrix would work
system(command =
         paste(path_to_tassel,
               "-importGuess", hmp,
               "-KinshipPlugin -method Centered_IBS -endPlugin -export ./two_step/kinship.txt -exportType SqrMatrix"))
kinship <- "./two_step/kinship.txt"

# run a mixed model to obtain residuals that will be used in the second step (GWAS)
system(command =
         paste(path_to_tassel,
               "-fork1 -r", pheno,
               "-fork2 -k", kinship,
               "-combine3 -input1 -input2 -MLMPlugin -endPlugin -export ./two_step/MLM"))

# TASSEL saves the residuals for each trait in a different file.
# The first N files are the residuals for N traits. The order they appear
# in the file is the order they will be saved. In this example, 
# MLM1 and MLM2 are the residuals for traits 1 and 2, respectively.
# For the next step, they need to be saved in a single file.
residuals1 <- fread(here("./two_step/MLM1.txt"), data.table = F, skip = 2)
residuals2 <- fread(here("./two_step/MLM2.txt"), data.table = F, skip = 2)
mt_residuals <- merge(residuals1, residuals2, sort = F)

# TASSEL requires phenotypic files to be saved with a specific head
# The commands below will save the multivariate residuals with the 
# required heads to run with the Multivariate Stepwise model below
write("<Phenotype>", "./two_step/mt_residuals.txt", sep = "\t")
fwrite(data.frame("taxa", "data", "data"), "./two_step/mt_residuals.txt", sep = "\t", append = T)
fwrite(mt_residuals, "./two_step/mt_residuals.txt", sep = "\t", append = T, row.names = F, quote = F)

# calling the multivariate stepwise plugin with an enter limit pvalue of 0.001.
# The code below runs the same model as in the manuscript. The only difference
# is that we are using residuals of a model fitted with a kinship matrix instead
# of the original phenotypes.
system(command =
         paste(path_to_tassel,
               "-fork1 -importGuess", hmp,
               "-fork2 -r ./two_step/mt_residuals.txt",
               "-combine3 -input1 -input2 -intersect",
               "-MultivariateStepwisePlugin  -usePerm false -enterLimit 0.001 -exitLimit 0.002", 
               "-saveToFile true -savePath ./two_step/two_step",
               "-parallel true -threads", ncores, "-endPlugin -export -runfork1 -runfork2"), intern = F)
