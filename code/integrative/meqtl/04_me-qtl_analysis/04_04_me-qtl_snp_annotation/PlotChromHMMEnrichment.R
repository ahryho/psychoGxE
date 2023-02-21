PlotChromHMMEnrichment <- function(enrichment.rslt, plt.title, plt.txt.size = 20){
  
  enrichment.rslt$state <- factor(enrichment.rslt$state, levels = enrichment.rslt$state)
  enrichment.rslt$state <- sub("[0-9]*_", "", enrichment.rslt$state)
  
  enrichment.rslt["is_sign"] <- ifelse(enrichment.rslt$p_val_emp < 0.05, "Significant", "Non-significant")
  
  ggplot(enrichment.rslt, aes(x = state, y = odds.ratio, fill = is_sign)) + 
    geom_bar( stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 1, color = "red") + 
    labs(x = "",
         y = "Odds ratio", 
         title = plt.title) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = plt.txt.size, colour = "black"),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 0.75 * plt.txt.size),
          axis.title = element_text(size = plt.txt.size),
          axis.text.x = element_text(size = 0.75 * plt.txt.size, angle = 90, colour = "black", hjust = 0.95, vjust = 0.2),
          axis.text.y = element_text(size = 0.75 * plt.txt.size, angle = 0, colour = "black")) +
    scale_fill_manual("", values = c( "#999999", "#009E73")) 
}
