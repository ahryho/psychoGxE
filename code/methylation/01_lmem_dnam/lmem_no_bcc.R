require(lme4)
require(foreign)

library(parallel)
library(foreach)
library(doParallel)

# LMEM with age, sex, BMI, and MDD status

# 1. Load data

args            <- commandArgs(T)
beta.mtrx.fn    <- as.character(args[1])
pheno.fn        <- as.character(args[2])
lmer.res.out.fn <- as.character(args[3])

beta.mtrx <- readRDS(beta.mtrx.fn) 
pheno     <- read.csv2(pheno.fn)
pheno     <- pheno[pheno$Include == 1, ]

pheno$Sample_ID  <- as.factor(pheno$Sample_ID)
pheno$Sex        <- as.factor(pheno$Sex)
pheno$Age        <- as.numeric(pheno$Age)
pheno$BMI        <- as.numeric(pheno$BMI)
pheno$Status     <- as.factor(pheno$Status)

# SVs:
pheno$DNAm_SV1  <- as.numeric(as.character(pheno$DNAm_SV1))
pheno$DNAm_SV2  <- as.numeric(as.character(pheno$DNAm_SV2), dec = ".")
pheno$DNAm_SV3  <- as.numeric(as.character(pheno$DNAm_SV3))

# Take only veh and dex sampl
samples.veh.ids <- pheno$DNAm_ID[pheno$Group == "veh"]
samples.dex.ids <- pheno$DNAm_ID[pheno$Group == "dex"]

# 2. Making sure about samples in pheno and and betas matrix in the same order

table(colnames(beta.mtrx) %in% pheno$DNAm_ID)

pheno   <- pheno[match(colnames(beta.mtrx), pheno$DNAm_ID ),]
all(pheno$DNAm_ID == colnames(beta.mtrx))

# 3. Build model

no.cores <- detectCores() - 1
cl <- makeCluster(no.cores)
registerDoParallel(cl)

res <- foreach(cpg =  1:nrow(beta.mtrx), .combine = rbind, .packages = 'lme4') %dopar% {

  sample <- pheno$Sample_ID 
  
  lmer.null  <- lmer(beta.mtrx[cpg, ] ~ pheno$Sex + pheno$Age + pheno$BMI + pheno$Status + 
                       (1|sample), REML = F)
  
  lmer.model <- lmer(beta.mtrx[cpg, ] ~ pheno$Sex + pheno$Age + pheno$BMI + pheno$Status + 
                       pheno$Group + (1|sample), REML = F)
  
  res.anova  <- anova(lmer.null, lmer.model)
  
  # Statistics
  mean.veh <- round(mean(beta.mtrx[cpg, samples.veh.ids]), 4)
  mean.dex <- round(mean(beta.mtrx[cpg, samples.dex.ids]), 4)
  sd.veh   <- round(sd(beta.mtrx[cpg, samples.veh.ids]), 4)
  sd.dex   <- round(sd(beta.mtrx[cpg, samples.dex.ids]), 4)
  var.veh  <- round(var(beta.mtrx[cpg, samples.veh.ids]), 4)
  var.dex  <- round(var(beta.mtrx[cpg, samples.dex.ids]), 4)
  
  fc       <- round(mean.dex - mean.veh, 4)# round(log2(mean.dex) - log2(mean.veh), 4)
  var.all  <- round(var(beta.mtrx[cpg, ]), 4)
  mad      <- round(mad(beta.mtrx[cpg, ]), 4)
  
  cbind(as.character(rownames(beta.mtrx)[cpg]), 
        signif(res.anova$"Pr(>Chisq)"[2], 4),
        round(res.anova$Chisq[2], 4),
        mean.veh, sd.veh, var.veh,
        mean.dex, sd.dex, var.dex,
        fc, var.all, mad)
}

stopImplicitCluster()

fdr <- p.adjust(res[, 2], method = "fdr")
res <- cbind(res, fdr)

res.colnames <- cbind("PROBE_ID", "PVAL", "CHI_SQ",
                      "MEAN_VEH", "SD_VEH", "VAR_VEH",
                      "MEAN_DEX", "SD_DEX", "VAR_DEX",
                      "FC", "VAR_ALL", "MAD", "FDR")
res          <- rbind(res.colnames, res)

write.table(res,
            file = lmer.res.out.fn, row.names = F, quote = F, sep = "\t", col.names = F, append = F)
