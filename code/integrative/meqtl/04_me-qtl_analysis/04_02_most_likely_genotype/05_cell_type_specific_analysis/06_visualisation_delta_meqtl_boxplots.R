setwd("~/bio/code/mpip/dex-stim-human-array/")

# Load functions

source("code/integrative/util.R")
source("code/integrative/meqtl/05_cell_type_enrichment/04_overlap_with_provencal_pnas_data/00_custom_functions.R")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)

## Load data
## 
dex.blood.types.snp.lst <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/dex_meqtl_snps_by_blood_cell_type_significant.rds") # dex == delta

# dex.blood.types.cpg.lst <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/dex_meqtl_cpgs_by_blood_cell_type_significant.rds") # dex == delta

delta.blood.types.meqtls.lst <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/delta_meqtls_by_blood_cell_type_significant.rds") 

delta.blood.types.meqtls.df <- delta.blood.types.meqtls.lst %>% bind_rows() %>% unique() %>% setDT()

# Extract meqtl for eqtm calculations
delta.meqtls.for.eqtms <- plyr::join(delta.blood.types.meqtls.df[, .(fdr = min(fdr)), by = CpG_ID], 
                                    delta.blood.types.meqtls.df[, .(CpG_ID, fdr, SNP)], match = "first")[, .(CpG_ID, SNP)]

write.table(delta.meqtls.for.eqtms, 
       "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/delta_meqtls_significant_for_eqtm_analysis.csv",
       row.names = F, quote = F, sep = "\t", col.names = T, append = F)


veh.blood.types.snp.lst <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/veh_meqtl_snps_by_blood_cell_type_significant.rds")

chromhmm.snp.state.lst.df <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/dex_meqtl_snps_by_chromhmm_state_blood_significant.rds")

epigenome.blood.map.df <- fread("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/epigenomes_blood.csv")

snp.df    <- fread("~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/snp_mtrx_backup/snp_mtrx.csv")

## DNAm beta matrices for delta meqtls
##
beta.veh.sub.mtrx <- GetMethylSubset(treatment = "veh", delta.meqtls.for.eqtms, is.out = F) # 653 x 197 
beta.dex.sub.mtrx <- GetMethylSubset(treatment = "dex", delta.meqtls.for.eqtms, is.out = F) # 653 x 197

## Load GEX matrices
## 
gex.dex.mtrx      <- LoadGEX("dex")
gex.veh.mtrx      <- LoadGEX("veh")

## Prepare data
## 
dex.blood.types.snps.df <- lapply(names(dex.blood.types.snp.lst), 
                                  function(blood.type){
                                    data.frame(Type = blood.type, SNP = dex.blood.types.snp.lst[[blood.type]])
                                  }) %>% bind_rows() %>% unique()

veh.blood.types.snps.df <- lapply(names(veh.blood.types.snp.lst), 
                                  function(blood.type){
                                    data.frame(Type = blood.type, SNP = veh.blood.types.snp.lst[[blood.type]])
                                  }) %>% bind_rows() %>% unique()

chromhmm.snps.df <- chromhmm.snp.state.lst.df %>% bind_rows()
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

## Function
## 
GetPlots <- function(blood.type, chromhmm.snps.df, dex.blood.types.snps.df, veh.blood.types.snps.df, chromhmm.state){
  
  if (blood.type %in% "Bmem+Bnv") {
    snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas == blood.type, ]$SNP), 
                    dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% c("Bmem", "Bnv"), "SNP"]),
                    baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in% c("Bmem", "Bnv"), "SNP"]))
  } else snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas %in% blood.type, ]$SNP), 
                         dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% blood.type, "SNP"]),
                         baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in%  blood.type, "SNP"]))
  
  intersect.dex.mesnp.chromhmm     <- intersect(snp.lst$chromHMM, snp.lst$dexDNAm)
  intersect.dex.veh.mesnp          <- intersect(snp.lst$baseDNAm, snp.lst$dexDNAm)
  intersect.dex.chromhmm.no.base   <- setdiff(intersect.dex.mesnp.chromhmm, intersect.dex.veh.mesnp)
  intersect.dex.chromhmm.with.base <- intersect(intersect.dex.mesnp.chromhmm, intersect.dex.veh.mesnp)
  
  intersect.veh.mesnp.chromhmm     <- intersect(snp.lst$chromHMM, snp.lst$baseDNAm)
  intersect.veh.chromhmm.no.dex    <- setdiff(intersect.veh.mesnp.chromhmm, intersect.dex.veh.mesnp)
  
  write.table(intersect.dex.chromhmm.no.base, 
              paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/for_phewas/",
                     "delta_mesnps_", blood.type, "_no_base.csv"),
              row.names = F, quote = F, sep = "\t", col.names = T, append = F)
  
  write.table(intersect.dex.chromhmm.with.base, 
              paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/for_phewas/",
                     "delta_mesnps_", blood.type, "_with_base.csv"),
              row.names = F, quote = F, sep = "\t", col.names = T, append = F)
  
  if (blood.type == "Bmem+Bnv") type <- c("Bmem", "Bnv") else type <- blood.type
  
  meqtl.id <- delta.blood.types.meqtls.df[SNP %in% intersect.dex.chromhmm.no.base,][Type %in% type][order(fdr, decreasing = F)] #[fdr == min(fdr)]
  plot.title <- paste0("Example of meQTL in ", blood.type, " in ", chromhmm.state, ", dex and chromhmm intersection\n")
  plt.no.base <- ProcessGetBoxPlot(beta.veh.sub.mtrx, beta.dex.sub.mtrx, snp.df, meqtl.id[1, ], fdr.thr = 0.05,  plot.title = plot.title)
  
  meqtl.id <- delta.blood.types.meqtls.df[SNP %in% intersect.dex.chromhmm.with.base,][Type %in% type][order(fdr, decreasing = F)] # [fdr == min(fdr)]
  plot.title <- paste0("Example of meQTL in ", blood.type, " in ", chromhmm.state, ", dex, baseline, chromhmm intersection\n")
  plt.with.base <- ProcessGetBoxPlot(beta.veh.sub.mtrx, beta.dex.sub.mtrx, snp.df, meqtl.id[1,], fdr.thr = 0.05,  plot.title = plot.title)
  
  grid.arrange(arrangeGrob(plt.no.base, plt.with.base,  nrow = 1)) 
}

## Plots
## 
chromhmm.state <- "all significant states"

blood.cell.types.of.interest <- c("Bmem+Bnv", "CD4mem", "CD4nv", "CD8mem", "Neu")
chromhmm.state.of.interest   <- c("2_TssAFlnk", "3_TxFlnk", "4_Tx", "7_Enh", "14_ReprPCWk", "6_EnhG", "10_TssBiv", "5_TxWk" )

blood.type <- blood.cell.types.of.interest[3]

plt.lst <- lapply(blood.cell.types.of.interest, GetPlots, 
                     chromhmm.snps.df = chromhmm.snps.df[state %in% chromhmm.state.of.interest, ], 
                     dex.blood.types.snps.df = dex.blood.types.snps.df, 
                     veh.blood.types.snps.df = veh.blood.types.snps.df,
                     chromhmm.state = chromhmm.state)

plt.lst

## eQTMs
## 
eqtms.fn <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/eqtm_cis_result_from_delta_meqtls_cell_type_fdr_005.csv"
eqtm.df <- fread(eqtms.fn)

selected.eqtm  <- eqtm.df[1,] 
selected.meqtl <- delta.blood.types.meqtls.df[CpG_ID == selected.eqtm$CpG_ID & SNP == selected.eqtm$SNP]
selected.eqtl  <- data.frame(ENSG_ID = selected.eqtm$ENSG_ID, SNP = selected.meqtl$SNP) %>% unique()

ScatterPlotGEXvsDNAm(beta.dex.sub.mtrx, beta.veh.sub.mtrx,
                     gex.dex.mtrx, gex.veh.mtrx,
                     cpg.id = selected.eqtm$CpG_ID, ensg.id = selected.eqtm$ENSG_ID)


ProcessGetBoxPlot(beta.veh.sub.mtrx, beta.dex.sub.mtrx, snp.df, selected.meqtl[1], fdr.thr = 0.05, 
                  plot.title = "")

ProcessGetBoxPlot(gex.veh.mtrx, gex.dex.mtrx, snp.df, selected.eqtl[1,], fdr.thr = 0.05, 
                  plot.title = "")
