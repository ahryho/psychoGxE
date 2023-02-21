
setwd("~/bio/code/mpip/dex-stim-human-array/")

# Load functions

source("code/integrative/util.R")
source("code/integrative/meqtl/05_cell_type_enrichment/04_overlap_with_provencal_pnas_data/00_custom_functions.R")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)

## Function
GetVennDiagram <- function(blood.type, chromhmm.snps.df, dex.blood.types.snps.df, veh.blood.types.snps.df, chromhmm.state){
  
  if (blood.type == "Bmem+Bnv") {
    snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas == blood.type, ]$CpG_ID), 
                    dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% c("Bmem", "Bnv"), "CpG_ID"]),
                    baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in% c("Bmem", "Bnv"), "CpG_ID"]))
  } else snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas %in% blood.type, ]$CpG_ID), 
                         dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% blood.type, "CpG_ID"]),
                         baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in%  blood.type, "CpG_ID"]))
  
  intersect.dex.meqtls   <- intersect(snp.lst$chromHMM, snp.lst$dexDNAm)
  perc.olap.dex.chromhmm <- scales::percent(length(intersect.dex.meqtls) / length(snp.lst$chromHMM), accuracy = 0.1)
  perc.olap.dex.dnam     <- scales::percent(length(intersect.dex.meqtls) / length(snp.lst$dexDNAm), accuracy = 0.1)
  
  intersect.veh.meqtls   <- intersect(snp.lst$chromHMM, snp.lst$baseDNAm)
  perc.olap.veh.chromhmm <- scales::percent(length(intersect.veh.meqtls) / length(snp.lst$chromHMM), accuracy = 0.1)
  perc.olap.veh.dnam     <- scales::percent(length(intersect.veh.meqtls) / length(snp.lst$baseDNAm), accuracy = 0.01)
  
  
  plot.title = paste0(paste0("Number of intersecting CpGs in ", blood.type, " in ", chromhmm.state, "\n"))
  
  cbPalette <- c("darkgrey", "cadetblue", "skyblue")
  
  ggVennDiagram(snp.lst, 
                category.names = c(paste0(perc.olap.dex.chromhmm, " - ChromHMM - ", perc.olap.veh.chromhmm), 
                                   paste0("dex\n", perc.olap.dex.dnam),
                                   paste0("base\n", perc.olap.veh.dnam)), 
                set_size = 4, 
                label_alpha = 0.7,
                edge_size = 0,
                set_geom = "text",
                set_color = "black",
                label = "count") +
    theme(legend.position = "none", 
          legend.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 12, face="italic")) +
    labs(x = "", y = "",
         title = plot.title) +
    # scale_color_brewer(palette = "Paired")
    # scale_fill_manual(values = alpha(cbPalette, 0.4))
    scale_fill_gradient(low = alpha(cbPalette[2], 0.7), high = alpha(cbPalette[3], 0.15))
}

## Load data
## 
dex.blood.types.cpg.lst <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/dex_meqtl_cpgs_by_blood_cell_type_significant.rds")

veh.blood.types.cpg.lst <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/veh_meqtl_cpgs_by_blood_cell_type_significant.rds")

chromhmm.cpg.state.lst.df <- readRDS("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/dex_meqtl_cpgs_by_chromhmm_state_blood_significant.rds")

epigenome.blood.map.df <- fread("~/bio/code/mpip/dex-stim-human-array/data/annotation/chromHMM/epigenomes_blood.csv")

##

dex.blood.types.cpgs.df <- lapply(names(dex.blood.types.cpg.lst), 
                                  function(blood.type){
                                    data.frame(Type = blood.type, CpG_ID = dex.blood.types.cpg.lst[[blood.type]])
                                    }) %>% bind_rows() %>% unique()

veh.blood.types.cpgs.df <- lapply(names(veh.blood.types.cpg.lst), 
                                  function(blood.type){
                                    data.frame(Type = blood.type, CpG_ID = veh.blood.types.cpg.lst[[blood.type]])
                                  }) %>% bind_rows() %>% unique()

cell.types.of.interest <- c("Bmem+Bnv", "CD4mem", "CD4nv", "CD8mem", "Neu")

# TssA (1)
# 
chromhmm.state <- "TssA"

chromhmm.snps.df <- chromhmm.cpg.state.lst.df$`1_TssA`
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- chromhmm.snps.df$Salas %>% unique()
blood.cell.types.of.interest <- blood.cell.types.of.interest[blood.cell.types.of.interest %in% cell.types.of.interest]

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df, dex.blood.types.snps.df = dex.blood.types.cpgs.df,
                     veh.blood.types.snps.df = veh.blood.types.cpgs.df, chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]], ncol = 2))

# TssAFlnk (2)
# 
chromhmm.state <- "TssAFlnk"

chromhmm.snps.df <- chromhmm.cpg.state.lst.df$`2_TssAFlnk`
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- chromhmm.snps.df$Salas %>% unique()
blood.cell.types.of.interest <- blood.cell.types.of.interest[blood.cell.types.of.interest %in% cell.types.of.interest]

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df, dex.blood.types.snps.df = dex.blood.types.cpgs.df,
                     veh.blood.types.snps.df = veh.blood.types.cpgs.df, chromhmm.state = chromhmm.state)

ggvenn.lst

# TxWk (5)
# 
chromhmm.state <- "TxWk"

chromhmm.snps.df <- chromhmm.cpg.state.lst.df$`5_TxWk`
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- chromhmm.snps.df$Salas %>% unique()
blood.cell.types.of.interest <- blood.cell.types.of.interest[blood.cell.types.of.interest %in% cell.types.of.interest]

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df, dex.blood.types.snps.df = dex.blood.types.cpgs.df,
                     veh.blood.types.snps.df = veh.blood.types.cpgs.df, chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]], ncol = 2), 
             arrangeGrob(ggvenn.lst[[3]], ggvenn.lst[[4]], ncol = 2),
             nrow = 2) 

# EnhG (6)
# 
chromhmm.state <- "EnhG"

chromhmm.snps.df <- chromhmm.cpg.state.lst.df$`6_EnhG`
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- chromhmm.snps.df$Salas %>% unique()
blood.cell.types.of.interest <- blood.cell.types.of.interest[blood.cell.types.of.interest %in% cell.types.of.interest]

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df, dex.blood.types.snps.df = dex.blood.types.cpgs.df,
                     veh.blood.types.snps.df = veh.blood.types.cpgs.df, chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]], ggvenn.lst[[3]], ncol = 3), 
             arrangeGrob(ggvenn.lst[[4]], ggvenn.lst[[5]], ncol = 2),
             nrow = 2) 

# Enh (7)
# 
chromhmm.state <- "Enhancer"

chromhmm.snps.df <- chromhmm.cpg.state.lst.df$`7_Enh`
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- chromhmm.snps.df$Salas %>% unique()
blood.cell.types.of.interest <- blood.cell.types.of.interest[blood.cell.types.of.interest %in% cell.types.of.interest]

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df, dex.blood.types.snps.df = dex.blood.types.cpgs.df, 
                     veh.blood.types.snps.df = veh.blood.types.cpgs.df, chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]], ggvenn.lst[[3]], ncol = 3), 
             arrangeGrob(ggvenn.lst[[4]], ggvenn.lst[[5]], ncol = 2),
             nrow = 2) 

# ReprPCWk
# 
chromhmm.state <- "ReprPCWk"

chromhmm.snps.df <- chromhmm.snp.state.lst.df$`14_ReprPCWk`
chromhmm.snps.df <- left_join(chromhmm.snps.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- chromhmm.snps.df$Salas %>% unique()
blood.cell.types.of.interest <- blood.cell.types.of.interest[blood.cell.types.of.interest %in% cell.types.of.interest]

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df, dex.blood.types.snps.df = dex.blood.types.snps.df, 
                     veh.blood.types.snps.df = veh.blood.types.snps.df, chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]],  nrow = 1)) 

######
###### Overlap between two cell type-specific analysis
###### 

chromhmm.state <- "all significant states"

chromhmm.cpgs.df <- chromhmm.cpg.state.lst.df %>% bind_rows()
chromhmm.cpgs.df <- left_join(chromhmm.cpgs.df, epigenome.blood.map.df[, .(EID, Salas)], by = c("eid" = "EID"))

blood.cell.types.of.interest <- c("Bmem+Bnv", "CD4mem", "CD4nv", "CD8mem", "Neu")
chromhmm.state.of.interest   <- c("2_TssAFlnk", "3_TxFlnk", "4_Tx", "7_Enh", "14_ReprPCWk", "6_EnhG", "10_TssBiv", "5_TxWk" )

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df[state %in% chromhmm.state.of.interest, ], 
                     dex.blood.types.snps.df = dex.blood.types.snps.df, 
                     veh.blood.types.snps.df = veh.blood.types.snps.df,
                     chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]], ggvenn.lst[[3]], ncol = 3), 
             arrangeGrob(ggvenn.lst[[4]], ggvenn.lst[[5]], ncol = 2),
             nrow = 2) 


#######
####### dex meQTLs with decrease or gain of methylation after dex
#######

ind.meqtl.delta.df      <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/me-qtl_cis_indp_rw_delta_fdr_005.csv") 

# Extract CpGs

beta.veh.sub.mtrx <- GetMethylSubset(treatment = "veh", ind.meqtl.delta.df, is.out = F) # 2'769 x 196 
beta.dex.sub.mtrx <- GetMethylSubset(treatment = "dex", ind.meqtl.delta.df, is.out = F) # 2'769 x 196

rownames(beta.veh.sub.mtrx) <- beta.veh.sub.mtrx$CpG_ID

beta.diff.veh.dex        <- beta.veh.sub.mtrx[, -1] - beta.dex.sub.mtrx[, -1]
beta.diff.veh.dex.median <- apply(beta.diff.veh.dex, 1, function(cpg) median(as.numeric(cpg)))

ind.meqtl.delta.beta.val <- data.frame(CpG_ID = beta.veh.sub.mtrx$CpG_ID,
                                       diff_veh_dex = beta.diff.veh.dex.median)

ind.meqtl.delta.beta.val[["meth_level"]] <- ifelse(ind.meqtl.delta.beta.val$diff_veh_dex > 0, "decreased", "gain")

ind.meqtl.delta.beta.val <- left_join(ind.meqtl.delta.beta.val, ind.meqtl.delta.df[, .(CpG_ID, SNP)])
ind.meqtl.delta.beta.val <- left_join(ind.meqtl.delta.beta.val, dex.blood.types.snps.df)

ind.meqtl.delta.df[["meth_level"]] <- ifelse(ind.meqtl.delta.df$beta > 0, "decreased", "gain")
ind.meqtl.delta.df <- left_join(ind.meqtl.delta.df, dex.blood.types.snps.df)

PlotMethyLevelGain <- function(blood.type, ind.meqtl.delta.beta.val, text.size = 14){
  
  if (blood.type == "Bmem+Bnv") blood.type <- c("Bmem", "Bnv")
  
  ggplot(ind.meqtl.delta.beta.val[ind.meqtl.delta.beta.val$Type %in% blood.type,], aes(x = meth_level)) +
    geom_bar(aes(fill = meth_level)) +
    geom_text(stat = 'count', aes(label = ..count..), vjust = -1, size = text.size / 3) +
    scale_y_continuous(position = "right") + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          plot.title = element_text(size = text.size), 
          axis.title = element_text(size = text.size),
          axis.text.y = element_blank(),
          axis.text = element_text(size = text.size, colour = "black")) + 
    labs(title = "",
         x = "", y = "Nr. of ind GR-meQTLs") +
    scale_fill_manual(values = c("#DC0000B2", "lightgrey")) 
}

barplt.lst <- lapply(blood.cell.types.of.interest, PlotMethyLevelGain, ind.meqtl.delta.beta.val = ind.meqtl.delta.beta.val)

barplt.lst <- lapply(blood.cell.types.of.interest, PlotMethyLevelGain, ind.meqtl.delta.beta.val = ind.meqtl.delta.df)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], barplt.lst[[1]], ncol = 2))
grid.arrange(arrangeGrob(ggvenn.lst[[2]], barplt.lst[[2]], ncol = 2)) 
grid.arrange(arrangeGrob(ggvenn.lst[[3]], barplt.lst[[3]], ncol = 2)) 
grid.arrange(arrangeGrob(ggvenn.lst[[4]], barplt.lst[[4]], ncol = 2)) 
grid.arrange(arrangeGrob(ggvenn.lst[[5]], barplt.lst[[5]], ncol = 2))

## interpret the overlapping/non-overlapping delta meSNPs and baseline meSNPs
## 
meqtl.veh.df <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/me-qtl_cis_result_veh_fdr_005.csv",  col.names = c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "FDR")) 

meqtl.dex.df <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/me-qtl_cis_result_dex_fdr_005.csv", col.names = c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "FDR")) 

meqtl.veh.sub.df <- meqtl.veh.df[SNP %in% unique(veh.blood.types.snps.df$SNP), .(SNP, CpG_ID, beta, FDR)]
meqtl.dex.sub.df <- meqtl.dex.df[SNP %in% unique(veh.blood.types.snps.df$SNP), .(SNP, CpG_ID, beta, FDR)]

meqtl.veh.dex.sub.df <- inner_join(meqtl.veh.sub.df, meqtl.dex.sub.df, 
                                   suffix = c("_veh", "_dex"), by = c("CpG_ID", "SNP"))

meqtl.veh.dex.sub.df[["sign"]] <- sign(meqtl.veh.dex.sub.df$beta_veh * meqtl.veh.dex.sub.df$beta_dex)

table(meqtl.veh.dex.sub.df$sign)

table(x.df$sign)

x.df <- meqtl.veh.dex.sub.df
x.df[["y"]] <- "sign"
x.df$sign <- ifelse(x.df$sign == 1, "positive", "opposite")
ggplot(x.df, aes( y = y, fill = sign, colour = sign)) + 
  geom_bar(position = "stack", alpha = 0.5) + 
  stat_count(geom = "text",
             # stat = "count",
             aes(label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
             angle = 0,  vjust = 0.5, hjust = 1.1, colour = "black", cex = 4) +
  stat_count(geom = "text",
             aes(label = ..count..),
             angle = 0,  vjust = -1, hjust = 1.5, colour = "black", cex = 4) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(title = "Distribution of parallel and opposite GR-meQTLs", 
       y = "", 
       x = "")


snp.opposite <- meqtl.veh.dex.sub.df[sign == -1, SNP]


# Venn plot
# 

chromhmm.state <- "all significant states"

blood.cell.types.of.interest <- c("Bmem+Bnv", "CD4mem", "CD4nv", "CD8mem", "Neu")
chromhmm.state.of.interest   <- c("2_TssAFlnk", "3_TxFlnk", "4_Tx", "7_Enh", "14_ReprPCWk", "6_EnhG", "10_TssBiv", "5_TxWk" )

ggvenn.lst <- lapply(blood.cell.types.of.interest, GetVennDiagram, 
                     chromhmm.snps.df = chromhmm.snps.df[state %in% chromhmm.state.of.interest, ][SNP %in% snp.opposite], 
                     dex.blood.types.snps.df = dex.blood.types.snps.df[dex.blood.types.snps.df$SNP %in% snp.opposite,], 
                     veh.blood.types.snps.df = veh.blood.types.snps.df[veh.blood.types.snps.df$SNP %in% snp.opposite, ], 
                     chromhmm.state = chromhmm.state)

grid.arrange(arrangeGrob(ggvenn.lst[[1]], ggvenn.lst[[2]], ggvenn.lst[[3]], ncol = 3), 
             arrangeGrob(ggvenn.lst[[4]], ggvenn.lst[[5]], ncol = 2),
             nrow = 2) 
