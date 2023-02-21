library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

# Load functions

source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation/EnrichmentWithPermutation_FUN.R")
source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation/RunEnrichmentChromHMM.R")

# Set up parameters

out.dir.pre             <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/"
# out.dir.pre           <- "~/bio/code/mpip/dex-stim-human-array/"
gwas.cluster.out.dir    <- paste0(out.dir.pre, "data/public_data/GWAS/")
gwas.enrichment.res.dir <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/with_missingness/enrichment/snps/")

# Load meqtl SNP data, Genomic Ranges object

meqtl.veh.snp.gr   <- readRDS(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtl_veh_snps_with_maf_gr.rds"))
meqtl.delta.snp.gr <- readRDS(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtl_delta_snps_with_maf_gr.rds"))

# Run GWAS enrichment

# GWAS CD2
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/GWAS/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds"))
# gwas.gw.gr <- gwas.gr[elementMetadata(gwas.gr)[, "p_value"] <= 5e-8, ] 

nperm <- 1000

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_CD_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.cd <- gwas.enrich.perm.rslt

# GWAS ADHD
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "/PGC/pgc_ADHD_Demontis_2019_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_ADHD_Demontis_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.adhd <- gwas.enrich.perm.rslt

# GWAS ASD
gwas.gr    <- readRDS(paste0(out.dir.pre, "data/public_data/PGC/pgc_ASD_Grove_2019_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_ASD_Grove_2019_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.asd <- gwas.enrich.perm.rslt

# GWAS BPD
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "PGC/BIP_2021/BIP_2021_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_BP_2021_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.bpd <- gwas.enrich.perm.rslt

# GWAS MDD
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "/PGC/pgc_MDD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_MDD_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.mdd <- gwas.enrich.perm.rslt

# GWAS SCZ
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "PGC/SCZ_2022/SCZ_2022_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_SCZ_2022_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.scz <- gwas.enrich.perm.rslt

# GWAS BMI
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "BMI_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_BMI_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.bmi <- gwas.enrich.perm.rslt

# GWAS IBD
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "IBD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_IBD_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.ibd <- gwas.enrich.perm.rslt

# GWAS PD
gwas.gr    <- readRDS(paste0(gwas.cluster.out.dir, "PD_GR_p005.rds"))

gwas.enrich.perm.rslt <- EnrichmentWithPermutation(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   public = gwas.gr, 
                                                   nperm = nperm)
gwas.enrich.perm.rslt[["data"]] <- "GWAS_PD_p0_05"
gwas.enrich.perm.rslt[["n_perm"]] <- nperm
gwas.enrich.perm.rslt[["n_gwas_snps"]] <- length(gwas.gr)

gwas.pd <- gwas.enrich.perm.rslt

(rslt <- rbind(gwas.cd, gwas.adhd, gwas.asd, gwas.bpd, gwas.mdd, gwas.ibd, gwas.bmi, gwas.scz, gwas.pd))

# Save enrichment results 
write.csv2(rslt, 
           file = paste0(gwas.enrichment.res.dir, "meqtl_snps_GWAS_enrichment_perm_", nperm, ".csv"), 
           row.names = F, quote = F)