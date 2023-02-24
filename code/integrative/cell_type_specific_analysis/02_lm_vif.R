
########################################################################
## Title: Script for checking for multicollinarity between blood cell types using VIF
## Author: Anastasiia Hryhorzhevska
########################################################################

## Section: load packages
########################################################################
require(foreign)
require(data.table)
require(ggplot2)
require(dplyr)
library(car)

## Section: working directory
########################################################################
setwd("~/bio/code/mpip/dex-stim-human-array/")
source("~/bio/code/mpip/dex-stim-human-array/code/integrative/util.R")

## Section: load data
########################################################################

treatment <- "dex"

beta.mtrx.fn    <- paste0("data/integrative/matrixEQTL/methyl_beta_mtrx_", treatment, ".csv")
pheno.fn        <- "data/pheno/pheno_full_for_kimono.csv"
bcc.fn          <- "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv" # Salas BCCs
out.fn          <- paste0("output/data/integrative/cell_type_enrichment/dnam_cell_type_enrichment_pvals_", treatment, ".csv")

# Load and prepare methyl beta mtrx

beta.mtrx       <- fread(beta.mtrx.fn) 
cpg.ids         <- data.frame(beta.mtrx[, CpG_ID])
beta.mtrx       <- beta.mtrx[, -c("CpG_ID")]
beta.mtrx       <- as.matrix(beta.mtrx)

# Load and prepare pheno data

pheno           <- read.csv2(pheno.fn) #, na.strings = "#N/A") 
dnam.veh.ids    <- pheno[pheno$Include == 1 & pheno$Group == "veh", "DNAm_ID"]
pheno           <- pheno[pheno$Include == 1 & pheno$Group == treatment, ]

pheno$Age               <- as.numeric(pheno$Age)
pheno$BMI_D1            <- as.numeric(pheno$BMI_D1)
pheno$PC1               <- as.numeric(pheno$PC1)
pheno$PC2               <- as.numeric(pheno$PC2)  
pheno$DNAm_SmokingScore <- as.numeric(pheno$DNAm_SmokingScore)

# blood cell counts (BCC)

bcc.df <- read.csv2(bcc.fn)
bcc.df <- bcc.df[match(dnam.veh.ids, bcc.df$DNAm_ID,), ] # put in the same order bcc df as pheno df
cov.df <- cbind(bcc.df[, -1], pheno[, c("DNA_ID", "Sex", "Age", "BMI_D1", "Status", "DNAm_SmokingScore", "PC1", "PC2")])

## Section: check whether the samples in the same order for all dfs
########################################################################
samples.ids <- as.character(cov.df$DNA_ID)
table(colnames(beta.mtrx) %in% samples.ids)
all(samples.ids == colnames(beta.mtrx))

cov.df <- cov.df[, -which(colnames(cov.df) == "DNA_ID")]

## Section: Build Model
########################################################################
cpg <- 1 #cg26928153
  
## LM
## 
lm.full               <- lm(beta.mtrx[cpg, ] ~ ., data = cov.df)
lm.no.neu             <- lm(beta.mtrx[cpg, ] ~ ., data = cov.df[, -which(colnames(cov.df) %in% c("Neu"))])
lm.no.cd4mem          <- lm(beta.mtrx[cpg, ] ~ ., data = cov.df[, -which(colnames(cov.df) %in% c("CD4mem"))])
lm.no.cd4mem.cd8mem   <- lm(beta.mtrx[cpg, ] ~ ., data = cov.df[, -which(colnames(cov.df) %in% c("CD4mem", "CD8mem"))])
lm.no.cd4mem.cd8mem.b <- lm(beta.mtrx[cpg, ] ~ ., data = cov.df[, -which(colnames(cov.df) %in% c("CD4mem", "CD8mem", "Bmem", "Bnv"))])

## 
## VIFs
vif.no.neu          <- data.frame(VIF = vif(lm.no.neu), lm = "no Neu") %>% mutate(Variable = rownames(.))
vif.no.cd4mem       <- data.frame(VIF = vif(lm.no.cd4mem), lm = "no CD4mem") %>% mutate(Variable = rownames(.))
vif.no.cd4and8mem   <- data.frame(VIF = vif(lm.no.cd4mem.cd8mem), lm = "no CD4mem & CD8mem") %>% mutate(Variable = rownames(.))
vif.no.cd4and8mem.b <- data.frame(VIF = vif(lm.no.cd4mem.cd8mem.b), lm = "no CD4mem & CD8mem & B-cells") %>% mutate(Variable = rownames(.))
  
vif.all <- rbind(vif.no.neu, vif.no.cd4mem, vif.no.cd4and8mem, vif.no.cd4and8mem.b)

## Section: Visualise VIFs
########################################################################
vif.plt <- ggplot(vif.all[vif.all$Variable %in% colnames(cov.df[1:12]), ]) +
  geom_bar(stat = "identity", 
            aes(y = VIF, x = Variable, fill = lm, color = lm), 
            alpha = 0.95, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 5, linewidth = 0.5, colour = "black") + 
  geom_text(aes(0, 5, label = "y = 5",  hjust = -1, vjust = -1)) +
  theme_custom() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey", linewidth = 0.25),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.25)) +
  labs(title = "", # paste0("VIF values for CpG = ", cpg.ids[cpg, ]),
         x = "", y = "Variance Inflation Factors") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")

vif.plt

ggsave("~/bio/code/mpip/dex-stim-human-array/output/plots/manuscript/m4_f2_vif_plt.pdf",
       vif.plt, width = 180, height = 100, units = "mm", dpi = 600, scale = 1.5)
