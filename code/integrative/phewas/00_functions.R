#' Title
#'
#' @param df 
#' @param number_of_hits 
#' @param xlab 
#' @param ylab1 
#' @param ylab2 
#' @param cbPalette 
#'
#' @return
#' @export
#'
#' @examples
extract_snp_top_hits_plot = function(df, number_of_hits = 10, 
                                     xlab = "No. of Associations", # --number_of_hits-- is replaced by variable
                                     ylab1 = "Top --number_of_hits-- Associations",
                                     ylab2 = "Associated CpGs",
                                     cbPalette = c("grey", "#0072B2")) {
  
  ylab1 <- str_replace(ylab1, "--number_of_hits--", as.character(number_of_hits))
  
  plot_df = df %>%
    add_count(trait_category, name = "n") %>%
    group_by(SNP) %>%
    add_count(SNP, name = "n_phenotypes") %>%
    arrange(desc(n_phenotypes)) %>%
    select(SNP, cpg, exposure_super_group, n_phenotypes) %>%
    filter(duplicated(SNP) == FALSE) %>% 
    distinct() %>%
    ungroup() %>%
    slice_head(n = number_of_hits) %>%
    dplyr::mutate(SNP = factor(unique(SNP), levels = rev(.$SNP))) 
  
  snp_hits_plot <- ggplot(plot_df, aes(y = as.numeric(SNP), x = n_phenotypes)) +
    geom_point(aes(fill = exposure_super_group), col = "black", shape = 21, size = 3) +
    labs(y = ylab1, 
         x = xlab) +
    scale_y_continuous(
      breaks = 1:length(plot_df$SNP),
      minor_breaks = NULL,
      labels = levels(plot_df$SNP),
      sec.axis = sec_axis(~., breaks = 1:length(plot_df$SNP),
                          name = ylab2,
                          labels = plot_df$cpg[as.numeric(plot_df$SNP)])
    ) +
    theme_bw() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size=text.size),
          legend.key.width = unit(0.8, 'cm'),
          plot.background = element_blank()) +
    scale_fill_manual(values = alpha(cbPalette, 0.5))
  
  return(snp_hits_plot)
}

#' Title
#'
#' @param mr_res 
#' @param blood_type 
#' @param cbPalette 
#'
#' @return
#' @export
#'
#' @examples
get_snp_hits_per_category_plot <- function(mr_res, bl_type, cbPalette = c("grey", "#0072B2")){
  snp_hits_trait_category <- mr_res[trait_category != "", ][blood_type == bl_type] %>%
    group_by(trait_category) %>%
    group_split()
  
  # bind list with all data to plot snp hits for all traits
  snp_hits_trait_category <- c(list(mr.res.ao.fdr.sig[trait_category != "" ]), snp_hits_trait_category)
  
  # plot
  snp_hits_category_plots <- map(snp_hits_trait_category, extract_snp_top_hits_plot, 
                                 number_of_hits = 10, ylab1 = "", ylab2 = "")
  
  for(i in 1:length(snp_hits_category_plots)){
    snp_hits_category_plots[[i]] <- snp_hits_category_plots[[i]] + 
      ggtitle(c(paste("All traits:", bl_type, sep = " "), 
                sort(unique(mr_res[trait_category != ""]$trait_category)))[i]) +
      theme(plot.title = element_text(face = "bold"))
  }
  
  return(snp_hits_category_plots)
}
