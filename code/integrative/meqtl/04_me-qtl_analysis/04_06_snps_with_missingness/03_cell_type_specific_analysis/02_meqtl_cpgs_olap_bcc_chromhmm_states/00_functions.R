cbPalette <- c("grey", "#0072B2")

#' Title
#'
#' @param blood.type 
#' @param chromhmm.snps.df 
#' @param dex.blood.types.snps.df 
#' @param veh.blood.types.snps.df 
#' @param chromhmm.state 
#'
#' @return
#' @export
#'
#' @examples
GetVennDiagram <- function(blood.type, chromhmm.snps.df, dex.blood.types.snps.df, veh.blood.types.snps.df, chromhmm.state, plot.title = NULL){
  
  if (blood.type == "Bmem+Bnv") {
    snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas == blood.type, ]$CpG_ID), 
                    dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% c("Bmem"), ]$CpG_ID),
                    baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in% c("Bmem"), ]$CpG_ID))
  } else snp.lst <- list(chromHMM = unique(chromhmm.snps.df[chromhmm.snps.df$Salas %in% blood.type, ]$CpG_ID), 
                         dexDNAm = unique(dex.blood.types.snps.df[dex.blood.types.snps.df$Type %in% blood.type, ]$CpG_ID),
                         baseDNAm = unique(veh.blood.types.snps.df[veh.blood.types.snps.df$Type %in%  blood.type, ]$CpG_ID))
  
  # plot.title <- paste0("Number of intersecting CpGs in ", 
  #                      ifelse(blood.type == "Bmem+Bnv", "Bmem", blood.type), 
  #                      " in ", chromhmm.state, "\n")
  if(blood.type == "Bmem+Bnv") blood.type <- "Bmem"; plot.title <- blood.type
  
  cbPalette  <- c("#0072B2", "transparent", "grey")
  
  euler.tbl <- euler(snp.lst)
  
  # For visualisation purposes:
  if (blood.type != "Neu") euler.tbl$ellipses["dexDNAm", c("a", "b")] <- euler.tbl$ellipses["dexDNAm", c("a", "b")] * 1.3

  plot(euler.tbl,
       quantities = list(type = c("counts")),
       fills = alpha(cbPalette[c(2, 1, 3)], 0.5),
       labels = list(c("", "ChromHMM", ""), cex = 0.75),
       main = c(plot.title, cex = 1))
}

##
#' Title
#'
#' @param blood.type 
#' @param ind.meqtl.delta.beta.val 
#' @param text.size 
#' @param plot.title 
#'
#' @return
#' @export
#'
#' @examples
PlotMethyLevelGain <- function(blood.type, ind.meqtl.delta.beta.val, text.size = 8, plot.title = NULL){
  if (blood.type == "Bmem+Bnv") blood.type <- c("Bmem", "Bnv")
  if (is.null(plot.title)) plot.title <- blood.type
  
  ggplot(ind.meqtl.delta.beta.val[ind.meqtl.delta.beta.val$Type %in% blood.type,], aes(x = meth_level), fill = meth_level, color = meth_level) +
    geom_bar(aes(fill = meth_level), alpha = 0.5) +
    geom_text(stat = 'count', aes(label = ..count..), vjust = 1.5, size = text.size / 3) +
    # scale_y_continuous(position = "right") + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          plot.title = element_text(size = text.size), 
          axis.title = element_text(size = text.size),
          axis.text.y = element_blank(),
          axis.text = element_text(size = text.size, colour = "black"),
          axis.ticks.y = element_blank()) + 
    labs(title = plot.title,
         x = "", y = "Nr. GR-meQTLs") +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2")
}

#' Title
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
BoxPlotDNAm <- function(df, plot.size = 8, plot.title = NULL){
  ggplot(df, aes(x = Type, y = value, color = Treatment, fill = Treatment)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    theme(legend.position = "bottom", # c(0.9, 0.9), 
          legend.title = element_blank(), 
          legend.text = element_text(size = plot.size),
          panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          plot.title = element_text(size = plot.size), 
          axis.title = element_text(size = plot.size),
          axis.text = element_text(size = plot.size, colour = "black")) + 
    labs(title = plot.title,
         x = "", y = "DNAm beta values") +
    scale_color_manual(values = cbPalette) +
    scale_fill_manual(values = alpha(cbPalette, 0.25))
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
