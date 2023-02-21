library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(IRanges)
library(arules)
library(reshape2)
library(stringr)

source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation/02a_generate_gr_functions.R")

out.dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/"

# 1. Load data

meqtl.dex.fn   <- paste0(out.dir.pre, "clumped/me-qtl_cis_result_snps_with_na_dex_beta_fdr_005.csv")
meqtl.veh.fn   <- paste0(out.dir.pre, "clumped/me-qtl_cis_result_snps_with_na_veh_beta_fdr_005.csv")
meqtl.delta.fn <- paste0(out.dir.pre, "clumped/me-qtl_cis_result_snps_with_na_delta_beta_fdr_005.csv")

col.names <-  c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")

ind.meqtl.dex.df   <- fread(meqtl.dex.fn)
ind.meqtl.veh.df   <- fread(meqtl.veh.fn)
ind.meqtl.delta.df <- fread(meqtl.delta.fn)

# 2. Generate background snps (baseline + delta )

snps.dir             <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/snps/"
background.all.gr.fn <- paste0(snps.dir, "dex_geno_snps_with_maf_bins_gr.rds")


if (file.exists(background.all.gr.fn)) background.all.gr <- readRDS(background.all.gr.fn) else {
  
  # load MAFs
  mafs           <- fread(paste0(snps.dir, "dex_geno_maf_for_filtered_samples.afreq"), data.table = F)
  colnames(mafs) <- c("CHR", "SNP", "REF", "ALT", "MAF", "OBS_CT")
  
  # generate bins of MAF in 0.05 steps
  boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
  mafs$bin <- arules::discretize(mafs$MAF, method = "fixed", breaks = boundaries, labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))
  
  # load background SNPs
  snps <- fread(paste0(snps.dir, "Dex_genoData_SNPs.bim"), data.table = F, select = c(2, 4), col.names = c("SNP", "POS")) # snp chr pos
  
  # create gr and save the object
  background.all.gr <- GenerateGrangesObjecteQTL(input = snps, 
                                                 mdata = mafs, 
                                                 ofile = background.all.gr.fn)
}

# 3. Generate merging file for eQTL data from background SNPs 

background <- as.data.frame(background.all.gr)[, c("seqnames", "start", "snp_id", "bin")]
colnames(background) <- c("CHR", "POS", "SNP", "bin")

meqtl.delta.snp.gr <- GenerateGrangesObjecteQTL(input = ind.meqtl.delta.df,
                                                mdata = background,
                                                ofile = paste0(out.dir.pre, "annotation/snps/meqtl_delta_snps_with_maf_gr.rds"))

meqtl.veh.snp.gr   <- GenerateGrangesObjecteQTL(input = ind.meqtl.veh.df,
                                                mdata = background,
                                                ofile = paste0(out.dir.pre, "annotation/snps/meqtl_veh_snps_with_maf_gr.rds"))

meqtl.dex.snp.gr   <- GenerateGrangesObjecteQTL(input = ind.meqtl.dex.df,
                                                mdata = background,
                                                ofile = paste0(out.dir.pre, "annotation/snps/meqtl_dex_snps_with_maf_gr.rds"))

# 4. Generate GR objects for GWAS Summary Statistics

gwas.cluster.out.dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/public_data/GWAS/"

# 4.1. Cross Disorders

fn <- "/binder/common/public_data/PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt" 

gwas.cross.dis <- fread(fn, 
                        header = T, stringsAsFactors = F, 
                        select = c("CHROM", "POS", "ID", "PVAL"))
colnames(gwas.cross.dis) <- c("CHR", "POS", "SNP", "P")

gwas.gr <- GenerateGrangesObject(gwas.cross.dis, 
                                 ofile = paste0(gwas.cluster.out.dir, "PGC/Cross_Disorder2/pgc_cdg2_meta_no23andMe_oct2019_v2_GR_p005.rds"))

# SCZ 

fn <- paste0(gwas.cluster.out.dir, "PGC/SCZ_2022/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz")

if(!file.exists(fn)){
  system(paste0("mkdir -p ", gwas.cluster.out.dir, "PGC/SCZ_2022/"))
  url <- "https://figshare.com/ndownloader/files/34517807/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz"
  download.file(url, fn)
}

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("ID", "CHROM", "POS", "PVAL"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "PGC/SCZ_2022/SCZ_2022_GR_p005.rds"))

# ADHD
 
fn <- "/binder/common/public_data/PGC/ADHD/Demontis_2019/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "pgc_ADHD_Demontis_2019_GR_p005.rds"))

# ASD

fn <- "/binder/common/public_data/PGC/ASD/Grove_2019/iPSYCHPGC_ASD_Nov2017"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "pgc_ASD_Grove_2019_GR_p005.rds"))

# Bipolar Disorder, BIP

fn <- paste0(gwas.cluster.out.dir, "PGC/BIP_2021/pgc-bip2021-BDI.vcf.tsv.gz")

if(!file.exists(fn)){
  system(paste0("mkdir -p ", gwas.cluster.out.dir, "PGC/BIP_2021/"))
  url <- "https://figshare.com/ndownloader/files/26603690/pgc-bip2021-BDI.vcf.tsv.gz"
  download.file(url, fn)
}

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("ID", "#CHROM", "POS", "PVAL"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "PGC/BIP_2021/BIP_2021_GR_p005.rds"))

# MDD

library(biomaRt)
ensembl <- useEnsembl("snp", dataset = "hsapiens_snp", GRCh = "37")

fn <- "/binder/common/public_data/PGC/MDD_2019/PGC_UKB_depression_genome-wide_ACGT.txt"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("MarkerName", "P"))
gwas <- gwas[P <= 0.05]

gwas.mdd.biomart <- getBM(attributes = c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
              filters = "snp_filter", values = gwas$MarkerName, mart = ensembl, uniqueRows = TRUE)

gwas <- left_join(gwas, gwas.mdd.biomart, by = c("MarkerName" = "SNP")) %>% 
  dplyr::select(SNP, CHR, POS, P)

colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenomicRanges::GRanges(seqnames = gwas$CHR,
                                  ranges = IRanges::IRanges(start = as.numeric(as.character(gwas$POS)),
                                                            end = as.numeric(as.character(gwas$POS))),
                                  snp_id = gwas$SNP,
                                  p_value = gwas$P)  

gwas.gr <- gwas.gr[!duplicated(gwas.gr)] 
GenomeInfoDb::seqlevelsStyle(gwas.gr) <- "UCSC"

saveRDS(gwas.gr, file =  paste0(gwas.cluster.out.dir, "pgc_MDD_GR_p005.rds"))

# IBD

fn <- "/binder/common/public_data/IBD/EUR.IBD.gwas_info03_filtered.assoc"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "IBD_GR_p005.rds"))

# PD

fn <- "/binder/common/public_data/PD/daner_panic_3GDES_col1_12"

gwas <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "CHR", "BP", "P"))
colnames(gwas) <- c("SNP", "CHR", "POS", "P")

gwas.gr <- GenerateGrangesObject(gwas, 
                                 ofile = paste0(gwas.cluster.out.dir, "PD_GR_p005.rds"))

# BMI

fn <- "/binder/common/public_data/BMI/Locke_et_al/SNP_gwas_mc_merge_nogc.tbl.uniq"

gwas.locke <- fread(fn, header = T, stringsAsFactors = F, select = c("SNP", "p"))

fn <- "/binder/common/public_data/BMI/ENGAGE/ENGAGE1000G_BMI.txt"

gwas <- fread(fn) #, header = T, stringsAsFactors = F, select = c("SNP", "chromosome", "position"))
gwas <- gwas %>% dplyr::select(SNP, CHR = chromosome, POS = position)

gwas.merged <- dplyr::left_join(gwas.locke, gwas)

colnames(gwas.merged) <- c("SNP", "P", "CHR", "POS")

gwas.gr <- GenerateGrangesObject(gwas.merged, 
                                 ofile = paste0(gwas.cluster.out.dir, "BMI_GR_p005.rds"))
