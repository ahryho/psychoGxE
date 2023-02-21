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
    scale_fill_gradient(low = alpha(cbPalette[2], 0.7), high = alpha(cbPalette[3], 0.15))
}

##

GetPlots <- function(blood.type, chromhmm.snps.df, dex.blood.types.snps.df, veh.blood.types.snps.df, chromhmm.state, eqtm.df,
                     beta.dex.sub.mtrx, beta.veh.sub.mtrx, bcc.df){
  
  if (blood.type %in% "Bmem+Bnv") {
    snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas == blood.type, ]$CpG_ID), 
                    dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% c("Bmem", "Bnv"), "CpG_ID"]),
                    baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in% c("Bmem", "Bnv"), "CpG_ID"]))
  } else snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas %in% blood.type, ]$CpG_ID), 
                         dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% blood.type, "CpG_ID"]),
                         baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in%  blood.type, "CpG_ID"]))
  
  intersect.dex.mesnp.chromhmm     <- intersect(snp.lst$chromHMM, snp.lst$dexDNAm)
  intersect.dex.veh.mesnp          <- intersect(snp.lst$baseDNAm, snp.lst$dexDNAm)
  intersect.dex.chromhmm.no.base   <- setdiff(intersect.dex.mesnp.chromhmm, intersect.dex.veh.mesnp)
  intersect.dex.chromhmm.with.base <- intersect(intersect.dex.mesnp.chromhmm, intersect.dex.veh.mesnp)
  
  intersect.veh.mesnp.chromhmm     <- intersect(snp.lst$chromHMM, snp.lst$baseDNAm)
  intersect.veh.chromhmm.no.dex    <- setdiff(intersect.veh.mesnp.chromhmm, intersect.dex.veh.mesnp)
  
  write.table(intersect.dex.chromhmm.no.base, 
              paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/for_phewas/cpg_level/",
                     "delta_mesnps_", blood.type, "_no_base.csv"),
              row.names = F, quote = F, sep = "\t", col.names = T, append = F)
  
  write.table(intersect.dex.chromhmm.with.base, 
              paste0("~/bio/code/mpip/dex-stim-human-array/output/data/integrative/cell_type_enrichment/for_phewas/cpg_level/",
                     "delta_mesnps_", blood.type, "_with_base.csv"),
              row.names = F, quote = F, sep = "\t", col.names = T, append = F)
  
  if (blood.type == "Bmem+Bnv") type <- c("Bmem", "Bnv") else type <- blood.type
  
  meqtl.id <- delta.blood.types.meqtls.df[CpG_ID %in% intersect.dex.chromhmm.no.base,][Type %in% type][order(fdr, decreasing = F)] #[fdr == min(fdr)]
  eqtm.id <- eqtm.df[CpG_ID %in% meqtl.id$CpG_ID][SNP %in% meqtl.id$SNP][FDR == min(FDR)]
  
  plot.title <- paste0("Example of meQTL in ", blood.type, " in ", chromhmm.state, ", dex and chromhmm intersection\n meQTL FDR = ", signif(meqtl.id[CpG_ID == eqtm.id$CpG_ID[1]][SNP == eqtm.id$SNP[1], fdr], 4))
  plt.no.base <- ProcessGetBoxPlot(beta.veh.sub.mtrx, beta.dex.sub.mtrx, snp.df, eqtm.id[1, ], fdr.thr = 0.05,  plot.title = plot.title)
  
  plot.title       <- paste0("Example of eQTM in ", blood.type, " in ", chromhmm.state, ", dex and chromhmm intersection\n eQTM FDR = ", signif(eqtm.id$FDR, 4))
  plt.eqtl.no.base <- ProcessGetBoxPlot(gex.veh.mtrx, gex.dex.mtrx, snp.df, eqtm.id[1, .(ENSG_ID, SNP)], fdr.thr = 0.05, 
                                        plot.title = plot.title)
  
  scatter.plt.no.base <- ScatterPlotGEXvsDNAm(beta.dex.sub.mtrx, beta.veh.sub.mtrx,
                                              gex.dex.mtrx, gex.veh.mtrx,
                                              cpg.id = eqtm.id$CpG_ID[1], ensg.id = eqtm.id$ENSG_ID[1])
  
  plt.bcc.dnam.no.base <- PlotBCCvsDNAm(beta.dex.sub.mtrx, beta.veh.sub.mtrx, bcc.df, cpg.id = eqtm.id$CpG_ID[1])
  
  # with base
  # 
  meqtl.id <- delta.blood.types.meqtls.df[CpG_ID %in% intersect.dex.chromhmm.with.base,][Type %in% type][order(fdr, decreasing = F)] # [fdr == min(fdr)]
  eqtm.id <- eqtm.df[CpG_ID %in% meqtl.id$CpG_ID][SNP %in% meqtl.id$SNP][FDR == min(FDR)]
  
  plot.title <- paste0("Example of meQTL in ", blood.type, " in ", chromhmm.state, ", dex, baseline, chromhmm intersection\n meQTL FDR = ", signif(meqtl.id[CpG_ID == eqtm.id$CpG_ID[1]][SNP == eqtm.id$SNP[1], fdr], 4))
  plt.with.base <- ProcessGetBoxPlot(beta.veh.sub.mtrx, beta.dex.sub.mtrx, snp.df, eqtm.id[1,], fdr.thr = 0.05,  plot.title = plot.title)
  
  plot.title <- paste0("Example of eQTM in ", blood.type, " in ", chromhmm.state, ", dex, baseline, chromhmm intersection\n eQTM FDR = ", signif(eqtm.id$FDR, 4))
  plt.eqtl.with.base <- ProcessGetBoxPlot(gex.veh.mtrx, gex.dex.mtrx, snp.df, eqtm.id[1, .(ENSG_ID, SNP)], fdr.thr = 0.05, 
                                          plot.title = plot.title)
  
  scatter.plt.with.base <- ScatterPlotGEXvsDNAm(beta.dex.sub.mtrx, beta.veh.sub.mtrx,
                                                gex.dex.mtrx, gex.veh.mtrx,
                                                cpg.id = eqtm.id$CpG_ID[1], ensg.id = eqtm.id$ENSG_ID[1])
  
  plt.bcc.dnam.with.base <- PlotBCCvsDNAm(beta.dex.sub.mtrx, beta.veh.sub.mtrx, bcc.df, cpg.id = eqtm.id$CpG_ID[1])
  
  return(list(no_base = list(plt.meqtl = plt.no.base, plt.eqtl = plt.eqtl.no.base, 
                             plt.eqtm = scatter.plt.no.base, plt.bcc.dnam = plt.bcc.dnam.no.base),
              with_base = list(plt.meqtl = plt.with.base, plt.eqtl = plt.eqtl.with.base, 
                               plt.eqtm = scatter.plt.with.base, plt.bcc.dnam = plt.bcc.dnam.with.base)))
}

##
##
PlotMethyLevelGain <- function(blood.type, ind.meqtl.delta.beta.val, text.size = 12){
  if (blood.type == "Bmem+Bnv") blood.type <- c("Bmem", "Bnv")
  
  ggplot(ind.meqtl.delta.beta.val[ind.meqtl.delta.beta.val$Type %in% blood.type,], aes(x = meth_level)) +
    geom_bar(aes(fill = meth_level)) +
    geom_text(stat = 'count', aes(label = ..count..), vjust = 1.5, size = text.size / 3) +
    scale_y_continuous(position = "right") + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          plot.title = element_text(size = text.size), 
          axis.title = element_text(size = text.size),
          axis.text.y = element_blank(),
          axis.text = element_text(size = text.size, colour = "black")) + 
    labs(title = "",
         x = "", y = "Nr. GR-meQTLs") +
    scale_fill_manual(values = c("#DC0000B2", "lightgrey"))
}

PlotBCCvsDNAm <- function(meth.beta.dex.mtrx, meth.beta.veh.mtrx, bcc.df, cpg.id){
    
    cbPalette <- c("#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
    
    if (blood.type  == "Bmem+Bnv"){
      blood.type <- c("Bmem", "Bnv")
      bcc.sub.df <- bcc.df %>% dplyr::select("DNA_ID", all_of(blood.type)) %>% setDT()
      bcc.sub.df <- bcc.sub.df[, Bcells := Bmem + Bnv][,.(DNA_ID, Bcells)]
      blood.type <- "Bmem+Bnv"
    } else bcc.sub.df <- bcc.df %>% dplyr::select("DNA_ID", all_of(blood.type)) %>% setDT()
    
    # meth.beta.veh.mtrx <- beta.veh.sub.mtrx
    # meth.beta.dex.mtrx <- beta.dex.sub.mtrx
    
    plt.df <- data.frame(CpG = c(as.numeric(meth.beta.dex.mtrx[CpG_ID == cpg.id, -1]),
                                 as.numeric(meth.beta.veh.mtrx[CpG_ID == cpg.id, -1])), 
                         treatment = c(replicate(ncol(meth.beta.dex.mtrx) - 1, "dex"), 
                                       replicate(ncol(meth.beta.dex.mtrx) - 1, "veh")),
                         DNA_ID = colnames(meth.beta.dex.mtrx)[- 1])
    
    plt.df <- left_join(plt.df, bcc.sub.df)
    colnames(plt.df) <- c("CpG", "treatment", "DNA_ID", "BCC")
    
    plt.df %>%
      ggplot(aes(x = CpG, 
                 y = BCC, 
                 color = treatment)) +
      geom_point() +
      geom_smooth(method = "lm") +
      xlab(cpg.id) + ylab(blood.type) +
      theme( 
        panel.background = element_blank(),
        plot.title = element_text(size = 10),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 0, hjust = 0.5), 
        legend.position = "bottom", 
        legend.title = element_blank()) +
      scale_color_manual(values = cbPalette)
}
