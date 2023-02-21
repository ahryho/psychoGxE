Sys.setenv(BIOMART_CACHE="~/bio/caches/biomaRt")
Sys.setenv(ANNOTATION_HUB_CACHE="~/bio/caches/AnnotationHub")

# Function to load libraries

LoadPackages <- function(pkg.list){
  new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}

LoadPackagesBio <- function(pkg.list){
  new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg, dependencies = TRUE)
  suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}

theme_custom <- function(plt.size = 12){
  theme( panel.background = element_blank(),
         plot.title = element_text(size = plt.size),
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_text(size = plt.size, color = "black"),
         axis.text.x = element_text(angle = 0, hjust = 0.5, size = plt.size, colour = "black"), 
         axis.text.y = element_text(size = plt.size, colour = "black"),
         legend.position = "none",
         legend.text = element_text(color = "black", size = plt.size)) 
}

LoadMethylBeta <- function(treatment){
  
  if (treatment == "delta")
    df.fn <- "/binder/mgp/workspace/2020_DexStim_Array_Human//dex-stim-human-array/output/data/integrative/matrixEQTL/dnam_residuals/methyl_beta_mtrx_delta.csv" else{
    dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/"
    df.fn    <- paste0(dir.pre, "methyl_beta_mtrx_", treatment ,".csv")
  }
  
  print("Loading DNAm data: \n")
  df       <- read_delim_arrow(df.fn, delim = ";")
  
  return(df)
}

LoadGEX <- function(treatment){
  dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/"
  df.fn    <- paste0(dir.pre, "gex_mtrx_", treatment ,".csv")

    print("Loading GEX data: \n")
    df       <- read_delim_arrow(df.fn, delim = ";")
  
  return(df)
}

LoadGenotype <- function(){
  dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/"
  df.fn    <- paste0(dir.pre, "snp_mtrx.csv")
  
  print("Loading Genotype data: \n")
  df       <- read_delim_arrow(df.fn, delim = ";")
  
  return(df)
}

LoadOmics <- function(treatment, is.methyl.df = T, is.gex.df = T, is.snp.df = T){
  
  if(is.methyl.df == T) 
    methyl.df <- LoadMethylBeta(treatment)
  
  if(is.gex.df == T) 
    gex.df <- LoadGEX(treatment)
  
  if(is.snp.df == T) 
    snp.df <- LoadGenotype()
  
  return(c(methyl.df, gex.df, snp.df))
}

LoadPheno <- function(treatment = ""){
  dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/pheno/"
  df.fn    <- paste0(dir.pre, "pheno_full_for_kimono.csv")
  
  print("Loading Phenotype data: \n")

  df      <- read.csv2(df.fn)
  
  if (treatment == "") 
    df <- df[df$Include == 1,] else 
      df <- df[df$Include == 1 & df$Group == treatment, ]
  
  df$Sample_ID         <- as.factor(df$Sample_ID)
  df$Status            <- as.factor(df$Status)
  df$Dex               <- as.factor(df$Dex)
  df$Age               <- as.numeric(df$Age)
  df$BMI_D1            <- as.numeric(df$BMI_D1)
  df$PC1               <- as.numeric(df$PC1)
  df$PC2               <- as.numeric(df$PC2)  
  df$DNAm_SV1          <- as.numeric(df$DNAm_SV1)
  df$DNAm_SV2          <- as.numeric(df$DNAm_SV2) 
  df$DNAm_SV3          <- as.numeric(df$DNAm_SV3) 
  df$DNAm_SmokingScore <- as.numeric(df$DNAm_SmokingScore)
  
  return(df)
}

LoadBio <- function(treatment){
  dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/data/integrative/matrixEQTL/"
  df.fn    <- paste0(dir.pre, "bio_mtrx_methyl_", treatment ,".csv")
  
  print("Loading GEX data: \n")
  df       <- fread(df.fn) 
  
  return(df)
}


LoadBCC <- function(is.dex = F){
  if(is.dex == T)
    df.fn    <- "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc_dex.csv" else
      df.fn    <- "output/data/methylation/dex_stim_array_human_epidish_salas_bcc_rpc.csv"
  
  print("Loading Salas BCCs data: \n")
  df       <- read.csv2(df.fn) 
  
  return(df)
}

GetVennPlt <- function(meqtl.df, eqtm.df, plot.title = NULL, cbPal.col = "#0072B2" ){
  meqtl.cpgs <- meqtl.df$CpG_ID %>% unique()
  eqtm.cpgs  <- eqtm.df$CpG_ID %>% unique()
  
  intersect.cpgs <- intersect(meqtl.cpgs, eqtm.cpgs)
  
  cpgs <- list(meqtl = meqtl.cpgs, 
               eqtm  = eqtm.cpgs)
  
  eqtm.neg.cor.df   <- eqtm.df[CpG_ID %in% intersect.cpgs][beta_eqtm < 0]
  eqtm.pos.cor.df   <- eqtm.df[CpG_ID %in% intersect.cpgs][beta_eqtm > 0]  
  
  perc.eqtm.neg.cor <- round(nrow(eqtm.neg.cor.df) / nrow(eqtm.df[CpG_ID %in% intersect.cpgs]) * 100, 1)
  perc.eqtm.pos.cor <- round(nrow(eqtm.pos.cor.df) / nrow(eqtm.df[CpG_ID %in% intersect.cpgs]) * 100, 1)
  
  perc.olap.mecpgs.with.ecpgs <- round(length(intersect.cpgs) / length(meqtl.cpgs) * 100, 1)
  perc.olap.ecpgs.with.mecpgs <- round(length(intersect.cpgs) / length(eqtm.cpgs) * 100, 1)
  
  if (is.null(plot.title))
    plot.title <- paste0(perc.olap.mecpgs.with.ecpgs, "% of meCpGs overlap with eQTM CpGs \n",
                         perc.eqtm.neg.cor, "% of eQTMs (in intersection with meQTLs) are negatively correlated \n",
                         perc.eqtm.pos.cor, "% of eQTMs (in intersection with meQTLs) are positively correlated \n")
  
  cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  p <- ggVennDiagram(cpgs, 
                category.names = c(paste0("meQTLs, ", perc.olap.mecpgs.with.ecpgs, "%"), 
                                   paste0("eQTMs, ", perc.olap.ecpgs.with.mecpgs, "%")), 
                label_alpha = 0.7,
                edge_size = 0,
                set_geom = "text",
                set_color = "black",
                label = "count") +
        theme(legend.position = "none", 
            legend.title = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(size = 10, face="italic")) +
      labs(x = "", y = "",
           title = plot.title) +
      scale_fill_gradient(low = "white", high = cbPal.col) 
  
  return(list(cpgs = intersect.cpgs, venn.plot = p))
}

GetFullmeQTLdf <- function(meqtl.df, snp.loc, cpg.loc, fdr.thr = 0.05){
  meqtl.fltr.df <- meqtl.df[fdr < fdr.thr, ] %>% unique()
  
  meqtl.full.df <- left_join(meqtl.fltr.df, snp.loc) %>% mutate(pos_snp = pos) %>% 
    dplyr::select(-c(pos))
  meqtl.full.df <- left_join(meqtl.full.df, cpg.loc[, .(CpG_ID, chr, pos)]) %>% 
    mutate(pos_cpg = pos) %>% 
    dplyr::select(-pos)
  meqtl.full.df <- meqtl.full.df[, dist := pos_snp - pos_cpg]  
  meqtl.full.df <- meqtl.full.df[, strand := ifelse(dist < 0, "-", "+")]
  meqtl.full.df <- meqtl.full.df[, abs_dist := abs(dist)]
  
  meqtl.full.df
}

GetFulleQTMdf <- function(df, cpg.loc = NULL, ensg.loc = NULL, fdr.thr = 0.05, treatment, is.dist  = F){
  colnames(df) <- c("CpG_ID", "ENSG_ID", "beta", "t-stat", "p-value", "fdr")
  
  df <- df[fdr < fdr.thr, ] %>% unique()
  
  # if (isTRUE(is.dist) | (cpg.loc != NULL & ensg.loc != NULL)){
  #   df <- left_join(df, ensg.loc)
  #   df <- left_join(df, cpg.loc[, .(CpG_ID, chr, pos_cpg)])
  # }
  
  df[["treatment"]] <- treatment
  df[["eQTM_ID"]]   <- paste(df$CpG_ID, df$ENSG_ID, sep = "-")
  
  # meqtl.full.df <- meqtl.full.df[, dist := pos_snp - pos_cpg]  
  # meqtl.full.df <- meqtl.full.df[, strand := ifelse(dist < 0, "-", "+")]
  # meqtl.full.df <- meqtl.full.df[, dist := abs(dist)]
  
  df[, .(eQTM_ID, CpG_ID, ENSG_ID, beta, `p-value`, fdr, treatment)] %>% setDT()
}

# function for Scatter Plot

GetScatterPlot2 <- function(df, selected.meqtl = NULL, fdr.thr = 0.05, plot.title = NULL, cbPalette = NULL, plt.size = 12){
  
  if (is.integer(cbPalette))
    cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  ggplot(df, aes(x = beta, y = -log10(fdr), shape = treatment, color = treatment)) +
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    geom_hline(yintercept = -log10(fdr.thr), colour = "#990000", linetype = "dashed") +
    geom_point(alpha = 1, size = 0.5) +
    scale_x_continuous(labels = scientific) +
    theme_custom() +
    labs( x = "Effect size (FC)", 
          title = plot.title) +
    scale_colour_manual(values = cbPalette)
}

GetScatterPlot3 <- function(df, selected.meqtl = NULL, fdr.thr = 0.05, plot.title = NULL){
  
  cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  ggplot(df, aes(x = methyl_change, y = -log10(fdr), shape = treatment, color = treatment)) +
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    geom_hline(yintercept = -log10(fdr.thr), colour = "#990000", linetype = "dashed") +
    geom_point(alpha = 1.5, size = 1.2) +
    scale_x_continuous(labels = scientific) +
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) +
    labs( x = "Effect size (FC)", 
          title = plot.title) +
    scale_colour_manual(values = cbPalette)
}

# Function for Boxplot

ToBed <- function(df, output.fn, is.save = T){
  bed <- df %>% dplyr::select(chrom = chr, chromStart = pos, chromEnd = pos, name = PROBE_ID)
  rownames(bed) <- bed$PROBE_ID
  
  if (is.save == T)
    write.table(bed,
                output.fn, 
                sep = " ", quote = F, row.names = F, col.names = F)
  bed
}

## Scatter plot and Boxplots for significant meQTLs -->
  
# func-getbetavalues
GetBetValuesDF <- function(methyl.beta.df, snp.df, selected.qtl, treatment){
  beta.values.df    <- methyl.beta.df[CpG_ID %in% selected.qtl$CpG_ID, -1]
  beta.values.df    <- rbind(beta.values.df, snp.df[SNP %in% selected.qtl$SNP, -1])
  beta.values.df    <- data.frame(t(beta.values.df)) %>% setDT()
  colnames(beta.values.df) <- c("CpG", "SNP")
  beta.values.df$SNP <- as.factor(beta.values.df$SNP)
  beta.values.df[SNP == 0, SNP := "AA"]; beta.values.df[SNP == 1, SNP := "AB"]; beta.values.df[SNP == 2, SNP := "BB"] 
  
  beta.values.df$treatment <- treatment
  beta.values.df
}

GetGEXValuesDF <- function(gex.df, snp.df, selected.qtl, treatment){
  values.df    <- gex.df[ENSG_ID %in% selected.qtl$ENSG_ID, -1]
  values.df    <- rbind(values.df, snp.df[SNP %in% selected.qtl$SNP, -1])
  values.df    <- data.frame(t(values.df)) %>% setDT()
  colnames(values.df) <- c("ENSG", "SNP")
  values.df$SNP <- as.factor(values.df$SNP)
  values.df[SNP == 0, SNP := "AA"]; values.df[SNP == 1, SNP := "AB"]; values.df[SNP == 2, SNP := "BB"] 
  
  values.df$treatment <- treatment
  values.df
}

# Scatter Plot
GetScatterPlot <- function(meqtl.all.full.df, selected.meqtl = NULL, fdr.thr = 0.05, plot.title = NULL){
  
  cbPalette <- c( "#0072B2", "#009E73", "#E69F00", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9", "#999999")
  
  ggplot(meqtl.all.full.df, aes(x = beta, y = -log10(fdr), shape = treatment, color = treatment)) +
    geom_vline(xintercept = 0, colour = "#990000", linetype = "dashed") + 
    geom_hline(yintercept = -log10(fdr.thr), colour = "#990000", linetype = "dashed") +
    geom_point(alpha = 1.5, size = 1.2) +
    geom_label_repel(data = selected.meqtl,
                     aes(x = beta,
                         y = -log10(fdr),
                         label = meQTL_ID),
                     fontface = 'bold',
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.5, "lines"),
                     segment.color = 'grey50',
                     nudge_x = 0.05, nudge_y = 10, 
                     size = 3) +
    scale_x_continuous(labels = scientific) +
    # scale_y_continuous(trans = trans_reverser('log10')) +
    # labs(title = " ", y = "", x = " ") + 
    theme(legend.position = "bottom", 
          legend.title = element_blank(),
          #  panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8)) +
    labs( x = "Effect size (FC), MatrixEQTL", 
          title = plot.title) +
    scale_colour_manual(values = cbPalette)
}

# BoxPlot
GetBoxPlot <- function(beta.values.df, selected.meqtl, fdr.thr = 0.05, plot.labels = c("AA", "AB", "BB"), plot.title = NULL){
  
  cbPalette <- c("grey", "#0072B2")
  
  beta.values.df %>%
    ggplot(aes(y = CpG, x = SNP, fill = treatment, color = treatment)) +
    geom_boxplot(width = 0.5) +
    # scale_fill_viridis(discrete = TRUE, alpha = 0.5) +
    scale_x_discrete(labels = plot.labels) +
    theme_custom() + 
    theme( legend.position = "bottom", 
           legend.title = element_blank(),
           legend.text = element_text(size = 12)) +
    labs(y = paste0(selected.meqtl$CpG_ID), 
         x = paste0("\n", selected.meqtl$SNP),
         title = plot.title) +
    scale_fill_manual(values = alpha(cbPalette, 0.25)) + 
    scale_color_manual(values = cbPalette)
}

# BoxPlot main
ProcessGetBoxPlot <- function(methyl.beta.veh.df, methyl.beta.dex.df, snp.df, selected.meqtl, 
                              fdr.thr = 0.05, plot.title = NULL){
  
  colnames(methyl.beta.veh.df)[1] <- "CpG_ID"
  colnames(methyl.beta.dex.df)[1] <- "CpG_ID"
  selected.meqtl                  <- selected.meqtl %>% dplyr::select(SNP, everything())
  colnames(selected.meqtl)[2]     <- "CpG_ID"
  
  beta.values.dex.df <- GetBetValuesDF(methyl.beta.dex.df, snp.df, selected.meqtl, "dex")
  beta.values.dex.df[["treatment"]] <- "GR-induced"
  beta.values.veh.df <- GetBetValuesDF(methyl.beta.veh.df, snp.df, selected.meqtl, "veh")
  beta.values.veh.df[["treatment"]] <- "Baseline"
  beta.values.df     <- rbind(beta.values.dex.df, beta.values.veh.df)
  
  snp.ind.cnt.df <- beta.values.dex.df %>% count(SNP) # beta.values.df[ , .(count = count(SNP)), by = .(SNP, treatment)]
  snp.ind.cnt.df[, label := paste0(SNP, "\n(n = ", n, ")")]
  
  GetBoxPlot(beta.values.df, selected.meqtl, snp.ind.cnt.df$label, fdr.thr = 0.05, plot.title)
}

GetManhattanPlot <- function(df, fdr.thr = 0.05, is.anno = T){
  
  plt.size = 12
  new.df <- df %>% 
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(CHR_LEN = max(pos)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(CHR_CUM_POS = cumsum(as.numeric(CHR_LEN)) - CHR_LEN) %>%
    dplyr::select(-CHR_LEN) %>%
    # Add this info to the initial dataset
    left_join(df, ., by = c("chr" = "chr")) %>%
    # Add a cumulative position of each SNP
    arrange(chr, pos) %>%
    mutate(BP_CUM_POS = pos + CHR_CUM_POS)
  
  axis.df <- new.df %>% 
    group_by(chr) %>% 
    summarize(center = (max(BP_CUM_POS) + min(BP_CUM_POS)) / 2)
  
  if (is.anno){
    top.10.min.fdr <- unique(new.df[, c("PROBE_ID", "UCSC_RefGene_Name", "BP_CUM_POS", "FDR")])
    top.10.min.fdr <- new.df[order(new.df$FDR)][1:15] 
    top.10.min.fdr[["GENE"]] <- apply(top.10.min.fdr, 1, function(x) strsplit(x[["UCSC_RefGene_Name"]], ";")[[1]][1])
    top.10.min.fdr$GENE <- ifelse(is.na(top.10.min.fdr$GENE), "", top.10.min.fdr$GENE)
  } else { top.10.min.fdr <- data.frame(matrix(nrow = 0, ncol = ncol(new.df) + 1)); colnames(top.10.min.fdr) <- c(colnames(new.df), "GENE")}
  
   cbPalette <- c( "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
    
   ggplot(new.df, aes(x = BP_CUM_POS, y = -log10(fdr))) +
    geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 0.5) +
    geom_hline(yintercept = -log10(fdr.thr), linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = rep(cbPalette, 2)) + 
    scale_x_continuous(label = axis.df$chr, breaks = axis.df$center, expand = c(0, 0)) +
    ggrepel::geom_text_repel(data = top.10.min.fdr,
                             aes(label = paste(PROBE_ID, GENE, sep = "\n"), x = BP_CUM_POS), 
                             max.overlaps = 15,
                             size =  plt.size / 4, color = "black") +
    labs( x = "Chromosomes",
          y = "-log10(fdr)") +
    theme( panel.background = element_blank(),
           plot.title = element_text(size = 10),
           # axis.line.x = element_blank(),
           # axis.line.y = element_blank(),
           panel.border = element_blank(),
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           axis.title = element_text(size = plt.size),
           axis.text.x = element_text(angle = 0, hjust = 0.5, size = plt.size, colour = "black"), 
           axis.text.y = element_text(size = plt.size, colour = "black"),
           legend.position = "none") 
}

ScatterPlotGEXvsDNAm <- function(meth.beta.dex.mtrx, meth.beta.veh.mtrx,
                                 gex.dex.mtrx, gex.veh.mtrx,
                                 cpg.id, ensg.id){
  
  cbPalette <-c("grey", "#0072B2")
  
  plt.df <- data.frame(CpG = c(as.numeric(meth.beta.dex.mtrx[CpG_ID == cpg.id, -1]),
                               as.numeric(meth.beta.veh.mtrx[CpG_ID == cpg.id, -1])), 
                       ENSG = as.numeric(gex.dex.mtrx[ENSG_ID == ensg.id, -1],
                                         gex.veh.mtrx[ENSG_ID == ensg.id, -1]), 
                       treatment = c(replicate(ncol(gex.dex.mtrx) - 1, "GR-induced"), replicate(ncol(gex.veh.mtrx) - 1, "baseline"))
  )
  
  plt.df %>%
    ggplot(aes(x = CpG, 
               y = ENSG, 
               color = treatment)) +
    geom_point() +
    geom_smooth(method = "lm") +
    xlab(cpg.id) + ylab(ensg.id) +
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

HistoPlotCntPerProbes <- function(df){
  
  num.ensg.per.cpg <- df %>% 
    group_by(CpG_ID) %>% tally()
  
  plot.gex.cnt <- ggplot(num.ensg.per.cpg, aes(n)) + 
    geom_bar(position = position_dodge()) +
    stat_count(geom = "text", 
               aes(label = comma(..count.., accuracy = 1L)),
               position = position_dodge(1),  vjust = -1, colour = "black", cex = 2) +
    scale_x_continuous(limits = c(0, max(num.ensg.per.cpg$n) + 5), 
                       breaks = seq(0, max(num.ensg.per.cpg$n) + 5, 5), 
                       expand = c(0, 0)) + 
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 8),
      axis.ticks.y = element_blank()) +
    labs(title = "Number of GEX probes per CpG", 
         x = "Number of GEX probes", 
         y = "")
  
  # Number of GEX probes per CpG
  
  num.cpg.per.gex <- venn.eqtm.df %>% 
    group_by(ENSG_ID) %>% tally()
  
  plot.cpg.cnt <- ggplot(num.cpg.per.gex, aes(n)) + 
    geom_bar(position = position_dodge()) +
    stat_count(geom = "text", 
               aes(label = comma(..count.., accuracy = 1L)),
               position = position_dodge(1),  vjust = -1, colour = "black", cex = 2) +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          axis.ticks.y = element_blank()) +
    labs(title = "Number of CpGs per GEX probe", 
         x = "Number of CpG sites", 
         y = "")
  
  return(list(gex = plot.gex.cnt, cpg  = plot.cpg.cnt))
}


# Take overlaps / non-overlaps

get_all_overlaps <- function(delta_df, dex_df, veh_df){
  intersect(intersect(veh_df, delta_df), dex_df)
}

get_veh_delta_overlaps <- function(delta_df, dex_df, veh_df){
  intersect(intersect(veh_df, delta_df), dex_df)
}

annotate_with_chipseeker <- function(df, chr, pos_start, pos_end, term = "SNP"){
  
  gr   <- GenomicRanges::GRanges(seqnames = paste0("chr", df[[chr]]),
                                 ranges = IRanges::IRanges(start = as.numeric(as.character(df[[pos_start]])),
                                                           end = as.numeric(as.character(df[[pos_end]]))))
  names(gr) <- df[[term]]
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  
  gr <- unique(gr)
  
  gr.anno <- annotatePeak(unique(gr), 
                          TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                          annoDb = "org.Hs.eg.db")
}

get_relative_freq_pos_across_chrom <- function(df, chr.length.df){
  freq.df        <- tabyl(df, chr, show_missing_levels = FALSE)

  freq.df                   <- left_join(freq.df, chr.length.df, by = "chr")
  freq.df$chr               <- factor(freq.df$chr,levels = chr.order)
  freq.df[["freq_chr"]]     <- freq.df$chr_size / sum(freq.df$chr_size)
  freq.df[["nr_cpg_chr"]]   <- freq.df$n / freq.df$freq_chr
  freq.df[["freq_cpg_chr"]] <- freq.df$nr_cpg_chr / sum(freq.df$nr_cpg_chr)
  
  return(freq.df)
}