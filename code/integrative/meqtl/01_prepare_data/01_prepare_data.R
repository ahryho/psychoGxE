library(dplyr)
library(data.table)
library(gdsfmt)
library(SNPRelate)

# Set up parameters

setwd("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/")

output.eqtm.pre <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/"
output.eqtm.pre <- "integrative/matrixEQTL/"

# src.pheno.data.pre <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/"
# src.snps.data.pre  <-"~/bio/datasets/snps/"
# src.pheno.data.pre <- "/binder/mgp/datasets/2020_DexStim_Array_Human/pheno/"

# Load data

db_gene_cg   <- fread("mapping/mapping_cpg_gene_ensg_full.csv")
db_ilmn_gene <- fread("mapping/mapping_ilmn_ensg_gene.csv")
db_gene_snp <- fread("~/bio/datasets/kimono/mapping/mapping_snp_gene_distance.csv")

snp.bim            <- fread("snps/final_imputed_qc_snps/filtered_196_samples/dex_geno_imputed_maf.bim")
pheno              <- fread("pheno/pheno_full_for_kimono.csv", na.strings = c('#N/A', "NA"), dec = ",") %>% setDT()
pheno              <- pheno[Include == 1]    
methyl.mtrx        <- readRDS("methylation/dex_methyl_beta_combat_mtrx.rds")
methyl.mval.mtrx   <- readRDS("methylation/dex_methyl_mval_combat_mtrx.rds")
gex.mtrx.veh       <- fread("integrative/matrixEQTL/gex_mtrx_veh.csv")
gex.mtrx.dex       <- fread("integrative/matrixEQTL/gex_mtrx_dex.csv")

# DNA_IDs to keep for further analysis

sample.ids.df <- data.frame(unique(pheno$DNA_ID), unique(pheno$DNA_ID))

fwrite(sample.ids.df, 
       "~/bio/code/mpip/dex-stim-human-array/data/snps/dna_ids_to_keep_for_analysis.txt",
       quote = F, row.names = F, col.names = F, sep = "\t")

# Prepare and save the SNP matrix

snps_gds <- snpgdsOpen("snps/gds/dex_geno_imputed.gds")
geno_obj <- snpgdsGetGeno(snps_gds, with.id = T)
snpgdsClose(snps_gds)

snp_mtrx   <- geno_obj$genotype
colnames(snp_mtrx)   <- geno_obj$snp.id
rownames(snp_mtrx)   <- geno_obj$sample.id

t_snp_mtrx <- t(snp_mtrx)
dim(t_snp_mtrx)
t_snp_mtrx <- cbind(SNP = colnames(snp_mtrx), t_snp_mtrx)
dim(t_snp_mtrx)
colnames(t_snp_mtrx)
t_snp_mtrx[1, 1:10] 

methyl.mtrx <- fread("integrative/matrixEQTL/methyl_beta_mtrx_veh.csv")

# snp.mtrx.na <- t_snp_mtrx %>% select(SNP, colnames(methyl.mtrx)[-1])

order.idx    <- c(0, match(colnames(methyl.mtrx)[-1], colnames(t_snp_mtrx)[-1])) + 1
snp.mtrx.na  <- t_snp_mtrx[, order.idx]

all(colnames(snp.mtrx.na)[-1] == colnames(methyl.mtrx)[-1])

fwrite(snp.mtrx.na, 
       "integrative/matrixEQTL/snp_mtrx_with_na.csv",
       quote = F, row.names = F, sep = ";")

snp.mtrx <- fread("integrative/matrixEQTL/snp_mtrx.csv")

# Prepare methylation data

# Extract only baseline samples:
veh.ids <- pheno[Dex == 0 & !is.na(DNAm_ID), .(DNA_ID, DNAm_ID)]
dex.ids <- pheno[Dex == 1 & !is.na(DNAm_ID), .(DNA_ID, DNAm_ID)]

# methyl.mtrx <- methyl.mval.mtrx

methyl.mtrx.veh <- methyl.mtrx[, colnames(methyl.mtrx) %in% veh.ids$DNA_ID]
methyl.mtrx.dex <- methyl.mtrx[, colnames(methyl.mtrx) %in% dex.ids$DNA_ID]

# Match colnames from methylation beta mtrx to values in "veh.ids" table
methyl.mtrx.veh           <- methyl.mtrx.veh[,match(veh.ids$DNAm_ID, colnames(methyl.mtrx.veh))]
colnames(methyl.mtrx.veh) <- veh.ids$DNA_ID
methyl.mtrx.veh           <- methyl.mtrx.veh %>% data.frame()
methyl.mtrx.veh["CpG_ID"] <- rownames(methyl.mtrx.veh)
methyl.mtrx.veh           <- methyl.mtrx.veh %>% dplyr::select(CpG_ID, everything())

methyl.mtrx.dex           <- methyl.mtrx.dex[,match(dex.ids$DNAm_ID, colnames(methyl.mtrx.dex))]
colnames(methyl.mtrx.dex) <- dex.ids$DNA_ID
methyl.mtrx.dex           <- methyl.mtrx.dex %>% data.frame()
methyl.mtrx.dex["CpG_ID"] <- rownames(methyl.mtrx.dex)
methyl.mtrx.dex           <- methyl.mtrx.dex %>% dplyr::select(CpG_ID, everything())

fwrite(methyl.mtrx.veh, 
       paste0(output.eqtm.pre, "methyl_mval_mtrx_veh.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(methyl.mtrx.dex, 
       paste0(output.eqtm.pre, "methyl_mval_mtrx_dex.csv"),
       quote = F, row.names = F, sep = ";")

# methyl.mtrx.dex <- fread(paste0(output.eqtm.pre, "methyl_beta_mtrx_dex.csv"), select = dex.ids$DNA_ID)
# methyl.mtrx.veh <- fread(paste0(output.eqtm.pre, "methyl_beta_mtrx_veh.csv"), select = veh.ids$DNA_ID)

# Calculate the differences between veh and dex

# Residuals

lmer.res.out.fn <- "/home/ahryhorzhevska/mpip/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/"
methyl.mtrx.veh <- fread(paste0(lmer.res.out.fn, "dnam_residuals_veh.csv"))
methyl.mtrx.dex <- fread(paste0(lmer.res.out.fn, "dnam_residuals_dex.csv"), select = colnames(methyl.mtrx.veh))

# methyl.mtrx.veh <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/dnam_residuals_veh.csv")
# methyl.mtrx.dex <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/dnam_residuals_dex.csv")

all(rownames(methyl.mtrx.veh) == rownames(methyl.mtrx.dex))
order.idx  <- match(colnames(methyl.mtrx.dex), colnames(methyl.mtrx.veh))

# methyl.mtrx.dex <- methyl.mtrx.dex[, colnames(methyl.mtrx.veh)]
all(colnames(methyl.mtrx.veh) == colnames(methyl.mtrx.dex))

methyl.mtrx.delta <- methyl.mtrx.veh[,-1] - methyl.mtrx.dex[,-1] 
methyl.mtrx.delta <- cbind(methyl.mtrx.veh[,1], methyl.mtrx.delta)

fwrite(methyl.mtrx.delta, 
       paste0(lmer.res.out.fn, "methyl_beta_mtrx_delta.csv"),
       # paste0(output.eqtm.pre, "methyl_beta_mtrx_delta.csv"),
       quote = F, row.names = F, sep = ";")

# Preapare GEX data
veh.ids <- pheno[Dex == 0 & !is.na(DNAm_ID), .(DNA_ID, RNA_ID)]
dex.ids <- pheno[Dex == 1 & !is.na(DNAm_ID), .(DNA_ID, RNA_ID)]

colnames(gex.mtrx.veh)[1] <- "RNA_ID"
colnames(gex.mtrx.dex)[1] <- "RNA_ID"

gex.mtrx.veh <- left_join(gex.mtrx.veh, veh.ids) %>% dplyr::select(-RNA_ID)
gex.mtrx.dex <- left_join(gex.mtrx.dex, dex.ids) %>% dplyr::select(-RNA_ID)

gex.mtrx.veh.t <- t(gex.mtrx.veh) %>% data.frame()
colnames(gex.mtrx.veh.t) <- gex.mtrx.veh$DNA_ID
gex.mtrx.veh.t <- gex.mtrx.veh.t[, match(veh.ids$DNA_ID, colnames(gex.mtrx.veh.t))]
gex.mtrx.veh.t["ENSG_ID"] <- rownames(gex.mtrx.veh.t)
gex.mtrx.veh.t <- gex.mtrx.veh.t %>% dplyr::select(ENSG_ID, everything())

gex.mtrx.dex.t <- t(gex.mtrx.dex) %>% data.frame()
colnames(gex.mtrx.dex.t) <- gex.mtrx.dex$DNA_ID
gex.mtrx.dex.t <- gex.mtrx.dex.t[,match(dex.ids$DNA_ID, colnames(gex.mtrx.dex.t))]
gex.mtrx.dex.t["ENSG_ID"] <- rownames(gex.mtrx.dex.t)
gex.mtrx.dex.t <- gex.mtrx.dex.t %>% dplyr::select(ENSG_ID, everything())

gex.mtrx.veh.t <- gex.mtrx.veh %>% dplyr::select(c("ENSG_ID", veh.ids$DNA_ID))
gex.mtrx.dex.t <- gex.mtrx.dex %>% dplyr::select(c("ENSG_ID", dex.ids$DNA_ID))

fwrite(gex.mtrx.veh.t, 
       paste0(output.eqtm.pre, "gex_mtrx_veh.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(gex.mtrx.dex.t, 
       paste0(output.eqtm.pre, "gex_mtrx_dex.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with SNP coordinates
snp.loc <- snp.bim %>% dplyr::select(V2, V1, V4) %>% unique()
colnames(snp.loc) <- c("SNP", "chr", "pos")

fwrite(snp.loc, 
       paste0(output.eqtm.pre, "snp_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with CpG coordinates

cpg.anno <- fread("~/bio/datasets/methylation/20_DMA/02_dmp/dmp_bcc_pcs_anno.csv")
colnames(cpg.anno)

cpg.loc           <- cpg.anno %>% dplyr::select(Name, chr, pos) %>% unique()
colnames(cpg.loc) <- c("CpG_ID", "chr", "pos")
cpg.loc$chr       <- gsub("[^0-9.]", "",  cpg.loc$chr)

order.idx  <- match(methyl.mtrx$CpG_ID, cpg.loc$CpG_ID)
cpg.loc    <- cpg.loc[order.idx, ]

all(cpg.loc$CpG_ID == methyl.mtrx$CpG_ID)

fwrite(cpg.loc, 
       paste0(output.eqtm.pre, "cpg_locations.csv"),
       quote = F, row.names = F, sep = ";")

cpg.loc <- fread(paste0(output.eqtm.pre, "cpg_locations.csv"))
cpg.loc[["pos_end"]] <- cpg.loc$pos

fwrite(cpg.loc, 
       paste0(output.eqtm.pre, "cpg_locations_meqtls.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare file with ENSG coordinates

ilmn.anno <- fread("~/bio/datasets/gene_expression/v3_v4_sharedContent_QC.txt")
ilmn.loc  <- ilmn.anno %>% dplyr::select(`X.PROBE_ID`, Chr, P_start, P_end) %>% unique()
colnames(ilmn.loc) <- c("Illumina_ID", "chr", "P_start", "P_end")

map.ilmn.ensg <- fread("~/bio/datasets/kimono/mapping/mapping_ilmn_ens.csv")
ensg.loc <- inner_join(ilmn.loc, map.ilmn.ensg) %>% dplyr::select(ENSG_ID = Ensemble_ID, chr, P_start, P_end) %>% unique() 

fwrite(ensg.loc, 
       paste0(output.eqtm.pre, "ensg_locations.csv"),
       quote = F, row.names = F, sep = ";")

# Prepare mapping file with snp ids, cpgs, ilmn ids, ennsg ids

map.cpg.ensg.gene <- fread("~/bio/datasets/kimono/mapping/mapping_cpg_gene_ensg_full.csv")
map.cpg.gene.ensg.ilmn <- inner_join(map.cpg.ensg.gene, map.ilmn.ensg)

fwrite(map.cpg.gene.ensg.ilmn, 
       paste0(output.eqtm.pre, "mapping_cpg_gene_ensg_ilmn_full.csv"),
       quote = F, row.names = F, sep = ";")

fwrite(map.cpg.gene.ensg.ilmn, 
       "~/bio/datasets/kimono/mapping/mapping_cpg_gene_ensg_ilmn_full.csv",
       quote = F, row.names = F, sep = ";")

# Prepare bio data

# pheno <- na.omit(pheno) # ?
pheno$DNAm_SV1 <- as.numeric(pheno$DNAm_SV1)
pheno$DNAm_SV2 <- as.numeric(pheno$DNAm_SV2)
pheno$DNAm_SV3 <- as.numeric(pheno$DNAm_SV3)

covariates <- colnames(pheno)

cov.list <- c("DNA_ID",
              "Sex", "Status", "Age", "BMI_D1", "DNAm_SmokingScore",
              # "V1", "V2", "V3",
              "DNAm_SV1", "DNAm_SV2", "DNAm_SV3", 
              "PC1", "PC2")

# VEH
bio.mtrx <- pheno[Dex == 0 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

all(colnames(methyl.mtrx.veh)[-1] == colnames(bio.mtrx.t)[-1])

order.idx  <- c(0, match(colnames(methyl.mtrx.veh)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_veh.csv"),
       quote = F, row.names = F, sep = ";")

# DEX
bio.mtrx <- pheno[Dex == 1 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

order.idx  <- c(0, match(colnames(methyl.mtrx.dex)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

all(colnames(methyl.mtrx.dex)[-1] == colnames(bio.mtrx.t)[-1])

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_dex.csv"),
       quote = F, row.names = F, sep = ";")

# bio.mtrx <- fread(paste0(output.eqtm.pre, "bio_mtrx_methyl.csv"))

# DELTA
cov.list <- c()

bio.mtrx.t <- pheno[Dex == 0 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

all(colnames(snp.mtrx.veh)[-1] == colnames(bio.mtrx.t)[-1])

order.idx  <- c(0, match(colnames(snp.mtrx.veh)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_delta.csv"),
       quote = F, row.names = F, sep = ";")


# Bio layer for eQTMs
# 

cov.list <- c("DNA_ID",
              "Sex", "Status", "Age", "BMI_D1", "DNAm_SmokingScore",
              "V1", "V2", "V3",
              "DNAm_SV1", "DNAm_SV2", "DNAm_SV3", 
              "PC1", "PC2")
 
# VEH

methyl.mtrx.veh <- fread("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_veh.csv")

bio.mtrx   <- pheno[Dex == 0 & !is.na(DNAm_ID), ..cov.list]
bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)

colnames(bio.mtrx.t)[1] <-  "Feature"

all(colnames(methyl.mtrx.veh)[-1] == colnames(bio.mtrx.t)[-1])

order.idx  <- c(0, match(colnames(methyl.mtrx.veh)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_gex_veh.csv"),
       quote = F, row.names = F, sep = ";")

# DEX

bio.mtrx   <- pheno[Dex == 1 & !is.na(DNAm_ID), ..cov.list]
bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

methyl.mtrx.dex <- fread("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/methyl_beta_mtrx_dex.csv")

order.idx  <- c(0, match(colnames(methyl.mtrx.dex)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

all(colnames(methyl.mtrx.dex)[-1] == colnames(bio.mtrx.t)[-1])

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_gex_dex.csv"),
       quote = F, row.names = F, sep = ";")


# DELTA
cov.list <- c("DNA_ID",
              "Sex", "Status", "Age", "BMI_D1", "DNAm_SmokingScore",
              "PC1", "PC2")

bio.mtrx.t <- pheno[Dex == 0 & !is.na(DNAm_ID), ..cov.list]

bio.mtrx.t <- dcast(melt(bio.mtrx.t, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

all(colnames(gex.mtrx.veh)[-1] == colnames(bio.mtrx.t)[-1])

order.idx  <- c(0, match(colnames(gex.mtrx.veh)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_gex_delta.csv"),
       quote = F, row.names = F, sep = ";")
# Bio layer with DNAm BCCs:

cov.list <- c("DNA_ID",
              "Sex", "Status", "Age", "BMI_D1", "DNAm_SmokingScore",
              "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran",
              "PC1", "PC2")

# VEH
bio.mtrx <- pheno[Dex == 0 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

all(colnames(methyl.mtrx.veh)[-1] == colnames(bio.mtrx.t)[-1])

order.idx  <- c(0, match(colnames(methyl.mtrx.veh)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_veh_dnam_bcc.csv"),
       quote = F, row.names = F, sep = ";")

# DEX
bio.mtrx <- pheno[Dex == 1 & !is.na(DNAm_ID), ] %>% dplyr::select(cov.list)

bio.mtrx.t <- dcast(melt(bio.mtrx, id.vars = "DNA_ID"), variable ~ DNA_ID)
colnames(bio.mtrx.t)[1] <-  "Feature"

order.idx  <- c(0, match(colnames(methyl.mtrx.dex)[-1], colnames(bio.mtrx.t)[-1])) + 1
bio.mtrx.t <- bio.mtrx.t[, ..order.idx]

all(colnames(methyl.mtrx.dex)[-1] == colnames(bio.mtrx.t)[-1])

fwrite(bio.mtrx.t, 
       paste0(output.eqtm.pre, "bio_mtrx_methyl_dex_dnam_bcc.csv"),
       quote = F, row.names = F, sep = ";")
