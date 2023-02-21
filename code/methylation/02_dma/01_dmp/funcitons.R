library(ggplot2)

TuneFcandPPlot <- function(df, pval.list, fc.list, plot.title = " "){
  est.df           <- data.frame(matrix(NaN, nrow = length(pval.list), ncol = length(fc.list)))
  colnames(est.df) <- fc.list
  rownames(est.df) <- pval.list
  
  for (p.thr in pval.list)
    for (delta.fc in fc.list){
      dmp.sign.df  <- df[abs(FC) > delta.fc & FDR < p.thr, ]
      est.df[as.character(p.thr), as.character(delta.fc)] <- nrow(dmp.sign.df)
    }
  
  est.df <- rownames_to_column(est.df, var = "p-value")
  est.df <- est.df %>% pivot_longer(cols = colnames(est.df)[-1], names_to = "FC", values_to = "Nr_DMPs")
  
  est.df$`p-value` <- as.numeric(est.df$`p-value`)
  est.df$FC        <- est.df$FC
  
  g.plot <- 
    ggplot(est.df, aes(`p-value`, Nr_DMPs, color = FC, group = FC)) + 
      geom_point() + geom_line() +
      theme( panel.background = element_blank(),
           plot.title = element_text(size = 10),
           axis.title = element_text(size = 10),
           axis.text.x = element_text(angle = 0, hjust = 1), 
           legend.position = "bottom") +
      labs(title = plot.title, y = "Number of DMPs", x = "FDR BH")
  
  return(list(plot = g.plot, data = est.df))
}

VolcanoPlot <- function(df, pval, fc, model = "SV", cbPalette){
  dmps.sign.df <- df[abs(FC) > fc & FDR < pval,]
  
  nr.hypermethyl <- nrow(dmps.sign.df[FC > fc & FDR < pval,])
  nr.hypomethyl  <- nrow(dmps.sign.df[FC < fc & FDR < pval,])
  
  df[, "MethylStatus"] <- ifelse(!(df$PROBE_ID %in% dmps.sign.df$PROBE_ID), 0, 
                                 ifelse(df$FC > fc, 1, 2)) %>% 
    as.factor()
  
  ggplot(df, aes(y = -log10(FDR), x = FC, color = MethylStatus)) +
    geom_point(size = 0.5) +
    scale_colour_manual(values = alpha(cbPalette, 0.5)) +
    geom_vline(xintercept = fc, colour = "black", linetype = "dotted") + 
    geom_vline(xintercept = -fc, colour = "black", linetype = "dotted") +
    geom_hline(yintercept = -log10(pval), colour = "black", linetype = "dotted") +
    theme( panel.background = element_blank(),
           plot.title = element_text(size = 12),
           axis.title = element_text(size = 12),
           axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
           axis.text.y = element_text(size = 12),
           legend.position = "none") +
    labs(title = "", # paste0("Volcano plot: ", model, " model with FDR <= ", pval, " and abs(FC) >= ", fc),
         x = "Fold Change", y = "-log10 FDR") +
    annotate(geom = "text", x = fc + 0.05,  y = -log10(min(df$FDR))/1.5, label = paste0(nr.hypermethyl, " \nhyper DMPs"), color = "black", fontface = "bold") +
    annotate(geom = "text", x = -fc - 0.05,  y = -log10(min(df$FDR))/1.5, label = paste0(nr.hypomethyl, " \nhypo DMPs"), color = "black", fontface = "bold") +
    xlim (-0.15, .15) +
    ylim(0, 18)
}
