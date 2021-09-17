A multivariate multi-locus stepwise approach for conducting GWAS on
correlated traits (Fernandes et al., 2021)
================

This repository contains all of the scripts and files we used to conduct
the simulation study presented here. All analyses were conducted on a
Unix system. In particular, `R` functions such as `system()` and
`mclapply()` may need to be modified to run on other systems.

### Load Packages

The library `here` is being used to automatically detect the working
directory (which should be the main folder of this repository).  
`simplePHENOTYPES` v1.3.0 was used in this simulation. In case of any
future changes that affect the results of this study, v1.3.0 can be
found here
(<https://github.com/samuelbfernandes/simplePHENOTYPES/releases/tag/v1.3.0>).
The `parallel` package is used for parallelization, and the function
`detectCores()` is used to determine the number of cores available for
that.

``` r
library(data.table)
library(here) 
library(fs)
library(simplePHENOTYPES)
library(SNPRelate)
library(tidyverse)
library(ggplot2)
library(LDheatmap)
library(flextable)
library(officer)
library(gdtools)
library(parallel)
ncores <- detectCores()
```

### Paths to External Software

All external software used for this simulation can be called from `R`.
To do that, the correct path to the executable should be provided below.
In the example below, the path to `TASSEL` is provided with additional
commands to increase the amount of memory available for the software.
Here it is set to 10GB. This value should be adjusted depending on the
amount of memory available.  
`TASSEL` is available here <https://tassel.bitbucket.io/>

`Plink` is available here <https://www.cog-genomics.org/plink/2.0/>

`GEMMA` is available here <https://anaconda.org/bioconda/gemma>

``` r
path_to_tassel <- "'/Applications/TASSEL\ 5/run_pipeline.pl' -Xms10G -Xmx10G"
path_to_plink <- "~/plink2"
path_to_gemma <- "/opt/anaconda3/bin/gemma"
```

### Create Directories

To save space, both marker data used in this study were compressed and
save in the folder `data`. This folder will be decompressed and
additional folders will be created to save the simulation results.

``` r
unzip("./data.zip")
dir_create(here("simulation"))
dir_create(here("results"))
```

### Load GEMMA

One of GWAS methods used in this study is the multivariate linear mixed
model (mvLMM) implemented in `GEMMA`. One of the things `GEMMA`requires
to run this model is a kinship matrix. We obtain a kinship matrix with
`GEMMA's` option `-gk 1`. The source file below contains a wrap that
calls `GEMMA` from `R` using the function `system()`. This file contains
two functions: `gemma()` and `relatedness()`.

``` r
source(here("./scripts/01_gemma.R"))
```

### Obtain Kinship From Marker Data

One of the file formats accepted by `GEMMA` is `Plink` BED files. For
simplicity, we used `TASSEL` to convert from HapMap to VCF and `Plink`
to convert from VCF to BED files. Next, `GEMMA` was used to obtain the
kinship matrix. Once the kinship matrix was obtained, we proceeded with
an eigen decomposition with `GEMMA` to speed up GWAS using mvLMM.  
All of these steps were didactically included here, but another (more
efficient) option would be to obtain a kinship matrix straight from
`TASSEL`.

``` r
#maize, sample size = 500
relatedness(path_to_tassel,
            path_to_plink,
            path_to_gemma,
            data = "./data/maize/ames_500_mac5_ld09_imputed")

#maize, sample size = 2815
relatedness(path_to_tassel,
            path_to_plink,
            path_to_gemma,
            data = "./data/maize/ames_2815_mac5_ld09_imputed")

#soybean, sample size = 500
relatedness(path_to_tassel,
            path_to_plink,
            path_to_gemma,
            data = "./data/soybean/soy_500_mac5_ld09")

#soybean, sample size = 2815
relatedness(path_to_tassel,
            path_to_plink,
            path_to_gemma,
            data = "./data/soybean/soy_2815_mac5_ld09")
```

### Obtain GDS Files

For detecting quantitative trait nucleotides (QTNs), we used an
LD-threshold of `r^2 = 0.2`. We used the function
`SNPRelate::snpgdsLDMat()` to calculate the LD between QTN and
significant SNPs from GWAS. To calculate LD, `snpgdsLDMat` uses the GDS
format, which can be exported from the
`simplePHENOTYPES::create_phenotypes()` function. All additional
parameters (e.g., add\_effect = 0.1) are only used to run
`create_phenotypes()`, they have no effect here.

``` r
#maize, sample size = 500
create_phenotypes(
  geno_file = "./data/maize/ames_500_mac5_ld09_imputed.hmp.txt",
  add_effect = 0.1,
  add_QTN_num = 1,
  rep = 1,
  h2 = 0.5,
  model = "A",
  output_dir = "./data/maize/gds",
  home_dir = here(), 
  quiet = T,
  verbose = F,
  out_geno  = "gds"
)
file_move("./data/maize/gds/ames_500_mac5_ld09_imputed.gds",
          "./data/maize/ames_500_mac5_ld09_imputed.gds")
unlink("./data/maize/gds", recursive = T)

#maize, sample size = 2815
create_phenotypes(
  geno_file = "./data/maize/ames_2815_mac5_ld09_imputed.hmp.txt",
  add_effect = 0.1,
  add_QTN_num = 1,
  rep = 1,
  h2 = 0.5,
  model = "A",
  output_dir = "./data/maize/gds",
  home_dir = here(), 
  quiet = T,
  verbose = F,
  out_geno  = "gds"
)
file_move("./data/maize/gds/ames_2815_mac5_ld09_imputed.gds",
          "./data/maize/ames_2815_mac5_ld09_imputed.gds")
unlink("./data/maize/gds", recursive = T)

#soybean, sample size = 500
create_phenotypes(
  geno_file = "./data/soybean/soy_500_mac5_ld09.hmp.txt",
  add_effect = 0.1,
  add_QTN_num = 1,
  rep = 1,
  h2 = 0.5,
  model = "A",
  output_dir = "./data/soybean/gds",
  home_dir = here(), 
  quiet = T,
  verbose = F,
  out_geno  = "gds"
)
file_move("./data/soybean/gds/soy_500_mac5_ld09.gds",
          "./data/soybean/soy_500_mac5_ld09.gds")
unlink("./data/soybean/gds", recursive = T)

#soybean, sample size = 2815
create_phenotypes(
  geno_file = "./data/soybean/soy_2815_mac5_ld09.hmp.txt",
  add_effect = 0.1,
  add_QTN_num = 1,
  rep = 1,
  h2 = 0.5,
  model = "A",
  output_dir = "./data/soybean/gds",
  home_dir = here(), 
  quiet = T,
  verbose = F,
  out_geno  = "gds"
)
file_move("./data/soybean/gds/soy_2815_mac5_ld09.gds",
          "./data/soybean/soy_2815_mac5_ld09.gds")
unlink("./data/soybean/gds", recursive = T)
```

### Run Simulations

This script runs all the simulations (see Table 1 on the manuscript)
conducted in this study. It generates all multivariate phenotypes and
all also exports the respective marker data without the SNPs selected to
be the QTNs.

``` r
source("./scripts/02_simulation.R")
```

### Run Multivariate Multi-locus Stepwise Model

The script below runs the multivariate multi-locus stepwise model in all
settings simulated in this study. This was the only instance where both
null and non-null settings were analyzed.

``` r
source("./scripts/03_mt_stepwise.R")
```

### Run Univariate Multi-locus Stepwise Model

The script below runs the univariate multi-locus stepwise model in all
non-null settings simulated in this study.

``` r
source("./scripts/04_st_stepwise.R")
```

### Run Multivariate Single Locus (GEMMA)

The script below runs the multivariate linear mixed model implemented in
`GEMMA` in all settings simulated in this study.

``` r
source("./scripts/05_mt_gemma.R")
```

### Process False Positive Detection Rate

The script below is used to process all GWAS results to obtain the false
positives rate from each setting.

``` r
source("./scripts/06_false_positive_rate.R")
```

### Process True Positive Detection Rate

The script below is used to process all GWAS results to obtain the true
positives rate from each setting.

``` r
source("./scripts/07_true_positive_rate.R")
```

### Generating Figures and Tables

The script below is used to generate all figures (except for figure 1
which is a trivial figure generated on excel) and the supplementary
table.

``` r
source("./scripts/08_graphs.R")
```

### Example of a Two-Step Approach Using Kinship Information

The multivariate multi-locus stepwise model does not allow random
effects. However, `TASSEL` provides all the necessary tools for running
a two-step approach. In the first step, we obtain residuals from a mixed
linear model including kinship information. In the second step, we use
these residuals as regular phenotypes in the multivariate multi-locus
stepwise model. The script below provides all the steps to run such an
analysis.

``` r
source("./scripts/09_example_two_step.R")
```
