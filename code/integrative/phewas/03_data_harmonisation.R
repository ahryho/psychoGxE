# ------------------------------------
# Data Harmonisation------------------
# ------------------------------------


# 1 Preparation-----------------------

## Load packages
library("tidyverse")
library("TwoSampleMR")
# remotes::install_github("mrcieu/TwoSampleMR")
library("ieugwasr")
# remotes::install_github("MRCIEU/MRInstruments")
library("MRInstruments")

## Load data
## 
## dir.pre     <-"/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/"
dir.pre       <- "~/bio/code/mpip/dex-stim-human-array/"
phewas.dir    <- paste0(dir.pre, "data/integrative/phewas/8_cell_type_groups/")
exposure.data <- readRDS(paste0(phewas.dir, "exposure_data_clumped.rds"))
ao            <- readRDS(paste0(dir.pre, "/data/public_data/PheWAS/IEU_List_Selected_Outcomes.rds"))

## Save unique exposure SNPs
unique.exposure.snps <- unique(exposure.data$SNP)

# 2 Extract outcome data--------------

outcome.data <- extract_outcome_data(
  snps = unique.exposure.snps,
  outcomes = ao$id,
  proxies = F
)

## Save data

saveRDS(outcome.data, 
        file = paste0(phewas.dir, "outcome_data_clumped.rds"))

# outcome.data <- readRDS(file = paste0(phewas.dir, "outcome_data_clumped.rds"))

# 3 Harmonisation---------------------

# exposure.data.no.base <- exposure.data[exposure.data$exposure %like% "no_base",]
# outcome.data.no.base  <- outcome.data[outcome.data$SNP %in% exposure.data.no.base$SNP,]

# exposure.data <- exposure.data.no.base # [1:8,]
# outcome.data  <- outcome.data.no.base # [outcome.data.no.base$SNP %in% exposure.data$SNP,]

data <- harmonise_data(
      exposure_dat = exposure.data,
      outcome_dat = outcome.data,
      action = 2
)

# 4 Bug fixing------------------------

## Some duplicate rows are present in the data. These are removed:
harmonised.data <- unique(data)

## The number of SNPs is added as a column
harmonised.data <- harmonised.data %>%
  add_count(exposure, outcome, name = "n_snp") %>%
  mutate(n_snp_1 = n_snp == 1)


# Looking at the instruments with >1SNP per instrument (i.e., SNP="rs9275314") shows that these outcomes have different effect sizes for the same SNP, so these are removed:
harmonised.data <- harmonised.data %>% 
  filter(n_snp_1 == TRUE) %>%
  select(-n_snp, n_snp_1)


# 5 Saving----------------------------

saveRDS(harmonised.data, 
        file = paste0(phewas.dir, "harmonised_data_clumped.rds"))

