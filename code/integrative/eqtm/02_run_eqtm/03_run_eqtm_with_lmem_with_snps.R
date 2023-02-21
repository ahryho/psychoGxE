library(data.table)
library(dplyr)

require(lme4)
require(foreign)

# 1. Load data

args         <- commandArgs(T)

treatment    <- as.character(args[1]) 
eqtm.in.pre  <- as.character(args[2]) 
eqtm.res.pre <- as.character(args[3]) 

mval <- "_beta"

eqtm.in.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/"
eqtm.res.pre <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/"

gex.layer.fn    <- paste0(eqtm.in.pre, "integrative/matrixEQTL/gex_mtrx_", treatment, ".csv")
methyl.layer.fn <- paste0(eqtm.in.pre, "integrative/matrixEQTL/methyl", mval, "_mtrx_", treatment, ".csv")
snp.layer.fn    <- paste0(eqtm.in.pre, "integrative/matrixEQTL/snp_mtrx.csv")

# Prepare mtrx

pheno.fn <- paste0(eqtm.in.pre, "pheno/pheno_full_for_kimono.csv")
pheno    <- read.csv2(pheno.fn)#, na.strings = "#N/A") 
pheno    <- pheno[pheno$Include == 1, ]

# Prepare GEX mtrx
# 
samples.ids <- pheno$RNA_ID[pheno$Group == treatment]

gex.layer  <- fread(gex.layer.fn)
colnames(gex.layer) <- c("ENSG_ID", samples.ids)

# Prepare methyl mtrx
# 
samples.ids <- pheno$DNAm_ID[pheno$Group == treatment]

dnam.layer  <- fread(methyl.layer.fn)
colnames(dnam.layer) <- c("CpG_ID", samples.ids)

# Prepare pheno
# 
pheno$Sample_ID  <- as.factor(pheno$Sample_ID)

pheno   <- pheno[match(colnames(dnam.layer)[-1], pheno$DNAm_ID ),]

# Load delta cell-type enriched meqtls
# 
delta.meqtls.for.eqtms <- fread(paste0(eqtm.res.pre, "/eqtms/delta_meqtls_significant_for_eqtm_analysis.csv"))

# Load genotype mtrx
# 
snp.layer   <- fread(snp.layer.fn)

# Load bio 
# 
bio.layer <-  pheno[, c("DNA_ID", "Sex", "Age", "BMI_D1", "Status", "DNAm_SmokingScore", "PC1", "PC2", 
                        "DNAm_SV1", "DNAm_SV2", "DNAm_SV3", "V1", "V2", "V3")]

# 2. Making sure about samples in pheno and and betas matrix in the same order

all(pheno$DNAm_ID == colnames(dnam.layer)[-1])
all(pheno$RNA_ID == colnames(gex.layer)[-1])
all(pheno$DNA_ID == colnames(snp.layer)[-1])
all(pheno$DNA_ID == bio.layer$DNA_ID)

# 3. Build model

res <- lapply(1:nrow(delta.meqtls.for.eqtms), function(meqtl){
  sample <- pheno$Sample_ID
  cpg    <- delta.meqtls.for.eqtms$CpG_ID[meqtl]
  snp    <- delta.meqtls.for.eqtms$SNP[meqtl]
  
  dnam <- t(as.matrix(dnam.layer[CpG_ID == cpg, -1]))
  geno <- t(as.matrix(snp.layer[SNP == snp, -1]))
  
  eqms.lst <- lapply(gex.layer$ENSG_ID, function(gene){
    gex  <- t(as.matrix(gex.layer[ENSG_ID == gene, -1]))
    
    df <- cbind(gex, dnam, geno, bio.layer[, -1])
    
    lm.model <- lm(df$gex ~ ., data = df)
    mdl.coef <- summary(lm.model)
    
    cbind(cpg, 
          gene, 
          p_value = signif(mdl.coef$coefficients["dnam","Pr(>|t|)"], 4),
          beta = signif(mdl.coef$coefficients["dnam","Estimate"], 4),
          snp)
  }) 
  
  data.frame(matrix(unlist(eqms.lst), nrow = length(eqms.lst), byrow=TRUE))
})

res <- bind_rows(res)
colnames(res) <- c("CpG_ID", "ENSG_ID", "PVAL", "FC", "SNP")

res[["FDR"]] <- p.adjust(res[, "PVAL"], method = "fdr")

eqtm.cis.result.fn   <- paste0(eqtm.res.pre, "eqtms/eqtm_cis_result_from_delta_meqtls_cell_type.csv")

write.table(res,
            file = eqtm.cis.result.fn, row.names = F, quote = F, sep = "\t", col.names = T, append = F)
