library(dplyr)
library(data.table)

# Set up parameters

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/")

output.eqtm.pre <- "integrative/matrixEQTL/"

# Load data

pheno         <- fread("pheno/pheno_full_for_kimono.csv", na.strings = c('#N/A', "NA"), dec = ",")[Include == 1] 
methyl.mtrx   <- fread("integrative/matrixEQTL/methyl_beta_mtrx_veh.csv")
snp.mtrx      <- fread("integrative/matrixEQTL/snp_mtrx_probabilities.csv")
snp.loc       <- fread("integrative/matrixEQTL/snp_locations.csv")

# Prepare and save the SNP matrix

order.idx  <- c(0, match(colnames(methyl.mtrx)[-1], colnames(snp.mtrx)[-1])) + 1
snp.mtrx   <- snp.mtrx[, order.idx]

all(colnames(snp.mtrx)[-1] == colnames(methyl.mtrx)[-1])

fwrite(snp.mtrx, 
       "integrative/matrixEQTL/snp_mtrx_probabilities.csv",
       quote = F, row.names = F, sep = ";")

snp.mtrx <- fread("integrative/matrixEQTL/snp_mtrx.csv")

# Prepare file with SNP coordinates

order.idx  <- match(snp.mtrx$SNP_ID, snp.loc$SNP)
snp.loc    <- snp.loc[order.idx, ]

all(snp.loc$SNP == snp.mtrx$SNP_ID)

fwrite(snp.loc, 
       paste0(output.eqtm.pre, "snp_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Check if snp and bio sample ids in the same order

bio.mtrx <- fread("integrative/matrixEQTL/bio_mtrx_methyl_dex.csv")
all(colnames(snp.mtrx)[-1] == colnames(bio.mtrx)[-1])

bio.mtrx <- fread("integrative/matrixEQTL/bio_mtrx_methyl_veh.csv")
all(colnames(snp.mtrx)[-1] == colnames(bio.mtrx)[-1])
