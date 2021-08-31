library(simplePHENOTYPES)#v1.3.0
library(here)
library(fs)
#Simulating traits for the null setting (validation)
cor_res <- matrix(c(1, 0.5, 0.5, 1), 2)
unzip("./data.zip")
dir_create(here("simulation"))
#---- Simulation null settings 500 ------
create_phenotypes(
  geno_file = here("data/maize/ames_500_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 3412, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_M",
  home_dir = here(), 
  cor_res = cor_res,
  verbose = F,
  quiet = T)

create_phenotypes(
  geno_file = here("data/soybean/soy_500_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 57687, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_S",
  home_dir = here(),
  cor_res = cor_res,
  verbose = F,
  quiet = T)

create_phenotypes(
  geno_file = here("data/maize/ames_500_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 2656, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_M",
  home_dir = here(), 
  verbose = F,
  quiet = T)

create_phenotypes(
  geno_file = here("data/soybean/soy_500_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 3477, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_S",
  home_dir = getwd(),
  verbose = F,
  quiet = T)

#---- Simulation null settings 2815 ----
create_phenotypes(
  geno_file = here("data/maize/ames_2815_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 3412, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_M_2815",
  home_dir = here(), 
  cor_res = cor_res,
  verbose = F,
  quiet = T)

create_phenotypes(
  geno_file = here("data/soybean/soy_2815_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 57687, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_S_2815",
  home_dir = here(),
  cor_res = cor_res,
  verbose = F,
  quiet = T)

create_phenotypes(
  geno_file = here("data/maize/ames_2815_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 2656, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_M_2815",
  home_dir = here(), 
  verbose = F,
  quiet = T)

create_phenotypes(
  geno_file = here("data/soybean/soy_2815_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0,0),
  seed = 3477, 
  rep = 100,
  ntraits = 2,
  h2 = c(0,0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_S_2815",
  home_dir = getwd(),
  verbose = F,
  quiet = T)

#------- Simulation of LD and Pleiotropy ----------
popsize <- c(500, 2815)
species <- c("soy", "maize")
arch <- c("LD", "pleiotropic")
effect <- c(0.5, 0.5)
grid <- expand.grid(popsize, species, arch, stringsAsFactors = F)
colnames(grid) <- c("popsize", "species", "arch")
grid <- split(grid, 1:nrow(grid))
h2 <- c(0.3, 0.8)
for(i in grid) {
  if( i$arch == "LD"){
    ld_min <- 0.5
    ld_max <- 0.9
  } else {
    ld_min <- NULL
    ld_max <- NULL
  }
  create_phenotypes(
    geno_file = here(ifelse(
      i$species == "soy",
      paste0("data/soybean/soy_", i$popsize, "_mac5_ld09.hmp.txt"),
      paste0("data/maize/ames_", i$popsize, "_mac5_ld09_imputed.hmp.txt")
    )),
    add_QTN_num = 3,
    add_effect = effect,
    rep = 100,
    ntraits = 2,
    h2 = h2,
    model = "A",
    ld_min = ld_min,
    ld_max = ld_max,
    type_of_ld = "direct",
    ld_method = "r",
    architecture = i$arch,
    output_format = "multi-file",
    output_dir = paste0("./simulation/",i$arch,"_",i$species,"_", i$popsize),
    home_dir = getwd(),
    quiet = T,
    remove_QTN = T,
    out_geno  = "BED",
    constraints = list(maf_above = 0.2, maf_below = 0.3)
  )
}

#------- Simulation of Partial Pleiotropy ----------
## Selecting QTNs for maize
# create_phenotypes(
#   geno_file = here("./data/maize/ames_500_mac5_ld09_imputed.hmp.txt"),
#   add_QTN_num = 6,
#   add_effect = 0.1,
#   rep = 1,
#   h2 = 1,
#   model = "A",
#   output_dir = paste0("./simulation/QTNs"),
#   home_dir = getwd(),
#   quiet = T,
#   seed = 999,
#   constraints = list(maf_above = 0.2, maf_below = 0.3)
# )
# SNPs_maize <- read.table("./simulation/QTNs/Additive_QTNs.txt", header = T)$snp
QTN_maize <-
  list(add = list(
    trait1 = c("S2_54360087", "S7_106651946", "S6_87674294"),
    trait2 = c("S2_54360087", "S7_106651946", "S4_127925619"),
    trait3 = c("S2_54360087", "S1_116015247", "S1_52638614")
  ))
dir_delete(here("./simulation/QTNs/"))

## Selecting QTNs for soybean
# create_phenotypes(
#   geno_file = here("./data/soybean/soy_500_mac5_ld09.hmp.txt"),
#   add_QTN_num = 6,
#   add_effect = 0.1,
#   rep = 1,
#   h2 = 1,
#   model = "A",
#   output_dir = paste0("./simulation/QTNs"),
#   home_dir = getwd(),
#   quiet = T,
#   seed = 888,
#   constraints = list(maf_above = 0.2, maf_below = 0.3)
# )
# SNPs_soy <- read.table("./simulation/QTNs/Additive_QTNs.txt", header = T)$snp
QTN_soy <-
  list(add = list(
    trait1 = c("ss715592677", "ss715632181", "ss715624441"),
    trait2 = c("ss715592677", "ss715632181", "ss715629188"),
    trait3 = c("ss715592677", "ss715611669", "ss715624933")
  ))
dir_delete(here("./simulation/QTNs/"))

popsize <- c(500, 2815)
species <- c("soy", "maize")
arch <- "partially"
effect <- c(0.5, 0.5, 0.5)
grid <- expand.grid(popsize, species, arch, stringsAsFactors = F)
colnames(grid) <- c("popsize", "species", "arch")
grid <- split(grid, 1:nrow(grid))
h2 <- c(0.3, 0.8, 0.8)
for(i in grid) {
  create_phenotypes(
    geno_file = here(ifelse(
      i$species == "soy",
      paste0("data/soybean/soy_", i$popsize, "_mac5_ld09.hmp.txt"),
      paste0("data/maize/ames_", i$popsize, "_mac5_ld09_imputed.hmp.txt")
    )),
    QTN_list = ifelse(i$species == "soy", QTN_soy, QTN_maize),
    ntraits = 3,
    add_effect = effect,
    rep = 100,
    h2 = h2,
    model = "A",
    architecture = i$arch,
    output_format = "multi-file",
    output_dir = paste0("./simulation/",i$arch,"_",i$species,"_", i$popsize),
    home_dir = getwd(),
    quiet = T,
    remove_QTN = T,
    out_geno  = "BED"
  )
}
