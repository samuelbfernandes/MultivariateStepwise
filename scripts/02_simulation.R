#-------------------------------------------------------------------------------
#         simulation of all the settings evaluated in this study
#-------------------------------------------------------------------------------

#----- Simulation of the null scenarios -----
# Correlation matrix used for the null scenarios
cor_res <- matrix(c(1, 0.5, 0.5, 1), 2)

#---- Simulation null settings 500 correlation = 0.5 ------
create_phenotypes(
  geno_file = here("data/maize/ames_500_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 3412,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_M",
  home_dir = here(),
  cor_res = cor_res,
  verbose = F,
  quiet = T
)

create_phenotypes(
  geno_file = here("data/soybean/soy_500_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 57687,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_S",
  home_dir = here(),
  cor_res = cor_res,
  verbose = F,
  quiet = T
)

#---- Simulation null settings 500 correlation = 0 ------
create_phenotypes(
  geno_file = here("data/maize/ames_500_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 2656,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_M",
  home_dir = here(),
  verbose = F,
  quiet = T
)

create_phenotypes(
  geno_file = here("data/soybean/soy_500_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 3477,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_S",
  home_dir = here(),
  verbose = F,
  quiet = T
)

#---- Simulation null settings 2815 correlation = 0.5 ------
create_phenotypes(
  geno_file = here("data/maize/ames_2815_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 3412,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_M_2815",
  home_dir = here(),
  cor_res = cor_res,
  verbose = F,
  quiet = T
)

create_phenotypes(
  geno_file = here("data/soybean/soy_2815_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 57687,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_cor_S_2815",
  home_dir = here(),
  cor_res = cor_res,
  verbose = F,
  quiet = T
)

#---- Simulation null settings 2815 correlation = 0 ------
create_phenotypes(
  geno_file = here("data/maize/ames_2815_mac5_ld09_imputed.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 2656,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_M_2815",
  home_dir = here(),
  verbose = F,
  quiet = T
)

create_phenotypes(
  geno_file = here("data/soybean/soy_2815_mac5_ld09.hmp.txt"),
  add_QTN_num = 1,
  add_effect = c(0, 0),
  seed = 3477,
  rep = 100,
  ntraits = 2,
  h2 = c(0, 0),
  model = "A",
  output_format = "multi-file",
  output_dir = "./simulation/null_S_2815",
  home_dir = here(),
  verbose = F,
  quiet = T
)

#------------------- Simulation of the non-null scenarios-----------------------
#------- Simulation of LD and Pleiotropy ----------
#parameters to be used in the simulation
popsize <- c(500, 2815)
species <- c("soy", "maize")
scenario <- c("LD", "pleiotropic")
allelic_effect <- c(0.5, 0.5)
grid <- expand.grid(popsize, species, scenario, stringsAsFactors = F)
colnames(grid) <- c("popsize", "species", "scenario")
grid <- split(grid, 1:nrow(grid))
h2 <- c(0.3, 0.8)
# Loop over all the different settings
for (i in grid) {
  # ld_min and ld_max are only used in the LD scenario
  # otherwise, they are set to NULL
  if (i$scenario == "LD") {
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
      paste0("data/maize/ames_", i$popsize,"_mac5_ld09_imputed.hmp.txt")
    )),
    add_QTN_num = 3,
    add_effect = allelic_effect,
    rep = 100,
    ntraits = 2,
    h2 = h2,
    model = "A",
    ld_min = ld_min,
    ld_max = ld_max,
    type_of_ld = "direct",
    ld_method = "r",
    architecture = i$scenario,
    output_format = "multi-file",
    output_dir = paste0("./simulation/", i$scenario, "_", i$species, "_", i$popsize),
    home_dir = here(),
    quiet = T,
    remove_QTN = T,
    out_geno  = "BED",
    constraints = list(maf_above = 0.2,
                       maf_below = 0.3)
  )
}

#------- Simulation of Partial Pleiotropy ----------
#------------------------------------------
#  A first round of simulation is conducted 
#  to select SNPs to be used as QTN in the
#  partial pleiotropy scenario. The actual
#  SNPs selected are used in the QTN_maize
#  and QTN_soy objects. None of the lines
#  commented below need to be run.
#------------------------------------------
#---- Selecting QTNs for maize ----
# create_phenotypes(
#   geno_file = here("./data/maize/ames_500_mac5_ld09_imputed.hmp.txt"),
#   add_QTN_num = 6,
#   add_effect = 0.1,
#   rep = 1,
#   h2 = 1,
#   model = "A",
#   output_dir = paste0("./simulation/QTNs"),
#   home_dir = here(),
#   quiet = T,
#   seed = 999,
#   constraints = list(maf_above = 0.2, maf_below = 0.3)
# )
# SNPs_maize <- read.table("./simulation/QTNs/Additive_QTNs.txt", header = T)$snp
# dir_delete(here("./simulation/QTNs/"))
#
#---- Selecting QTNs for soybean -----
# create_phenotypes(
#   geno_file = here("./data/soybean/soy_500_mac5_ld09.hmp.txt"),
#   add_QTN_num = 6,
#   add_effect = 0.1,
#   rep = 1,
#   h2 = 1,
#   model = "A",
#   output_dir = paste0("./simulation/QTNs"),
#   home_dir = here(),
#   quiet = T,
#   seed = 888,
#   constraints = list(maf_above = 0.2, maf_below = 0.3)
# )
# SNPs_soy <- read.table("./simulation/QTNs/Additive_QTNs.txt", header = T)$snp
# dir_delete(here("./simulation/QTNs/"))

#parameters to be used in the simulation
#SNPs pre-selected to be QTNs when simulating partial pleiotropy in maize
QTN_maize <-
  list(add = list(
    trait1 = c("S2_54360087", "S7_106651946", "S6_87674294"),
    trait2 = c("S2_54360087", "S7_106651946", "S4_127925619"),
    trait3 = c("S2_54360087", "S1_116015247", "S1_52638614")
  ))
#SNPs pre-selected to be QTNs when simulating partial pleiotropy in soybean
QTN_soy <-
  list(add = list(
    trait1 = c("ss715592677", "ss715632181", "ss715624441"),
    trait2 = c("ss715592677", "ss715632181", "ss715629188"),
    trait3 = c("ss715592677", "ss715611669", "ss715624933")
  ))
popsize <- c(500, 2815)
species <- c("soy", "maize")
scenario <- "partially"
allelic_effect <- c(0.5, 0.5, 0.5)
grid <- expand.grid(popsize, species, scenario, stringsAsFactors = F)
colnames(grid) <- c("popsize", "species", "scenario")
grid <- split(grid, 1:nrow(grid))
h2 <- c(0.3, 0.8, 0.8)
# Loop over all the different settings
for (i in grid) {
  create_phenotypes(
    geno_file = here(ifelse(
      i$species == "soy",
      paste0("data/soybean/soy_", i$popsize, "_mac5_ld09.hmp.txt"),
      paste0("data/maize/ames_", i$popsize, "_mac5_ld09_imputed.hmp.txt")
    )),
    QTN_list = ifelse(i$species == "soy", QTN_soy, QTN_maize),
    ntraits = 3,
    add_effect = allelic_effect,
    rep = 100,
    h2 = h2,
    model = "A",
    architecture = i$scenario,
    output_format = "multi-file",
    output_dir = paste0("./simulation/", i$scenario, "_", i$species, "_", i$popsize),
    home_dir = here(),
    quiet = T,
    remove_QTN = T,
    out_geno  = "BED"
  )
}
