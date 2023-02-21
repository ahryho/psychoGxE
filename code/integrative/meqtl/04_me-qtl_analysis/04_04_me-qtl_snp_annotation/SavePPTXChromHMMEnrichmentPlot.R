# tissue <- "blood"
# 
# enrichment.rslts.fn <- paste0(out.dir.pre, "meqtl_snps_chromHMM_", tissue, "/", list.files(paste0(out.dir.pre, "meqtl_snps_chromHMM_", tissue)))
# 
# pptx.template <- read_pptx("~/bio/code/mpip/dex-stim-human-array/output/plots/enrichment/meqtl_snp_chromHMM_enrichment.pptx")
# pptx.out      <- paste0("~/bio/code/mpip/dex-stim-human-array/output/plots/enrichment/meqtl_snp_chromHMM_", tissue, "_enrichment.pptx")
# 

library(officer)

#' Function to save enrichment plots to .pptx file
#'
#' @param enrichment.rslts.fn a list of file names which contain enrichment results. Each file represent a cell type from ChromHMM
#' @param tissue blood or brain
#' @param epigenome an epigenome mapping table
#' @param pptx.template an rpptx object
#' @param pptx.out a presentation name which will be saved
#' 
#' @export pptx.out a pptx object 
#'
#' @examples
#' 

SavePPTXChromHMMEnrichmentPlot <- function(enrichment.rslts.fn, tissue, epigenome, pptx.template, pptx.out, plt.txt.size = 16){
  
  enrichment.all.df <- lapply(seq_along(enrichment.rslts.fn), function(i) {
    eid   <- sub(".csv", "", sub(".*_", "", enrichment.rslts.fn[i]))
    ename <- epigenome[EID == eid, NAME] 
    ename <- paste(eid, gsub("? from.*$", "", ename), sep = ": ")  
    # ename <- paste0(eid, ": ", ename)
    plot.title <- paste0("Histone mark enrichment for the delta meSNPs for ", ename)
    
    df           <- read.csv2(enrichment.rslts.fn[i])
    colnames(df) <- c("n_overlap", "odds.ratio", "or_perm", "p_val", "p_val_perm", "p_val_emp", "state", "n_perm")
    ggp <- PlotChromHMMEnrichment(df, plt.title = plot.title)
    
    add_slide(pptx.template,
              layout = "Plot",
              master = "Basic") %>%
      ph_with(
        value = ggp,
        location = ph_location_label(ph_label = "Content Placeholder 2")
      )
    return(data.frame(df, ename = ename))
  }) %>%
    bind_rows()
  
  enrichment.all.df$state <- factor(enrichment.all.df$state, levels = unique(enrichment.all.df$state))
  enrichment.all.df$state <- sub("[0-9]*_", "", enrichment.all.df$state)
  
  # Heatmap plot for emperical p-values
  
  ggp.heatmap.pval <- 
    ggplot(enrichment.all.df, aes(state, ename, fill = -log10(p_val_emp))) + 
      geom_tile(colour = "white", size = 0.15) +
      theme(legend.position = "top",
            legend.title = element_text(size = 0.75 * plt.txt.size, colour = "black"),
            legend.text = element_text(size = 0.75 * plt.txt.size, colour = "black"),
            legend.key.width = unit(2, 'cm'),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(size = 0.75 * plt.txt.size),
            axis.title = element_blank(),
            axis.text.x = element_text(size = 0.75 * plt.txt.size, angle = 90, colour = "black", hjust = 0.95, vjust = 0.2),
            axis.text.y = element_text(size = 0.75 * plt.txt.size, angle = 0, colour = "black")) +
      scale_fill_viridis(name = "-log10(P-value)")
  
  add_slide(pptx.template,
            layout = "Plot",
            master = "Basic") %>%
    ph_with(
      value = ggp.heatmap.pval,
      location = ph_location_label(ph_label = "Content Placeholder 2")
    )
  
  # Heatmap plot for odd ratios
  
  enrichment.all.df <- enrichment.all.df[order(enrichment.all.df$ename), ]
  
  ggp.heatmap.or <- 
    ggplot(enrichment.all.df, aes(state, ename, fill = odds.ratio)) + 
      geom_tile(colour = "white", size = 0.15) +
      theme(legend.position = "top",
            legend.title = element_text(size = 0.75 * plt.txt.size, colour = "black"),
            legend.text = element_text(size = 0.75 * plt.txt.size, colour = "black"),
            legend.key.width = unit(2, 'cm'),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(size = 0.75 * plt.txt.size),
            axis.title = element_blank(),
            axis.text.x = element_text(size = 0.75 * plt.txt.size, angle = 90, colour = "black", hjust = 0.95, vjust = 0.2),
            axis.text.y = element_text(size = 0.75 * plt.txt.size, angle = 0, colour = "black")) +
      scale_fill_viridis(limits = c(floor(min(enrichment.all.df$odds.ratio)), round(max(enrichment.all.df$odds.ratio))), 
                         breaks = seq(floor(min(enrichment.all.df$odds.ratio)), round(max(enrichment.all.df$odds.ratio)), by = 0.5))
  
  # write.csv2(enrichment.all.df, 
  #            file = "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_blood_enrichment_delta_vs_veh_by_type.csv", 
  #            row.names = F, quote = F)
  
  add_slide(pptx.template,
            layout = "Plot",
            master = "Basic") %>%
    ph_with(
      value = ggp.heatmap.or,
      location = ph_location_label(ph_label = "Content Placeholder 2")
    )
  
  print(pptx.template, target = pptx.out)
  
  return(list(pval = ggp.heatmap.pval, or = ggp.heatmap.or))
}
