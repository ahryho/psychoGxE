EnrichmentWithPermutation <- function(own, background, public, nperm){
  
   # own <- meqtl.delta.snp.gr
   # public <-  chromhmm.blood.states[(elementMetadata(chromhmm.blood.states)[, "type"]) == "5_TxWk", ] # chromhmm.all.states
   # background <- meqtl.veh.snp.gr
  
  GenomeInfoDb::seqlevelsStyle(own) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(public) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(background) <- "UCSC"
  
  overlap     <- sum(overlapsAny(own, public))  # number of overlaps between the delta snps and public data (chromHMM_state or GWAS)
  non_overlap <- length(own) - overlap
  
  overlap_bkgr      <- sum(overlapsAny(background, public))  # number of overlaps between baseline snps and public data (chromHMM_state or GWAS)
  non_overlap_bkgr  <- length(background) - overlap_bkgr
  
  conf_mtrx        <- matrix(c(overlap, overlap_bkgr, non_overlap, non_overlap_bkgr), 2, 2, byrow = TRUE)
  fisher_test_rslt <- fisher.test(conf_mtrx)
  p_value <- fisher_test_rslt$p.value
  or      <- fisher_test_rslt$estimate
  
  # resampling with respect to the MAF bins 
  
  background_bins <- lapply(1:11, function(x) background[background$bin == x])  # all background SNPs in MAF bin 1 to 11
  own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # length of all 11 MAF bins
  
  resampling <- lapply(1:nperm, function(x){
    sample_overlap <- sum(sapply(1:11, function(y)
      # check that the background maf bin is not empty
      ifelse(length(background_bins[[y]]) != 0, 
             # take n random samples of the background (n = own bin length except background bin length is smaller, than take all background genes of this MAF bin) 
             sum(IRanges::overlapsAny(sample(background_bins[[y]], 
                                             ifelse(length(background_bins[[y]]) >= own_bin_lengths[y], own_bin_lengths[y], length(background_bins[[y]]))), 
                                      public)), 
             0)
    )
    )
    sample_non_overlap <- length(own) - sample_overlap
    conf_mtrx <- matrix(c(overlap, sample_overlap, non_overlap, sample_non_overlap), 2, 2, byrow = TRUE)
    fisher_test_rslt <- fisher.test(conf_mtrx)
    # c(p_value = fisher_test_rslt$p.value, fisher_test_rslt$estimate)
    c(p_value = fisher_test_rslt$p.value, fisher_test_rslt$estimate, sample_overlap = sample_overlap)
  }
  )

  p_value_emp <- (sum(unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] <= p_value) + 1) / (nperm + 1)
  p_value_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "p_value"] %>% mean() 
  or_permutat <- unlist(resampling)[attr(unlist(resampling), "names") == "odds ratio"] %>% mean() 
  
  p_value_karo <- ifelse((sum(unlist(resampling)[attr(unlist(resampling), "names") == "sample_overlap"] >= overlap) / nperm) == 0, 
                         paste("<=", 1/nperm, sep = " "), 
                         sum(unlist(resampling)[attr(unlist(resampling), "names") == "sample_overlap"] >= overlap) / nperm)
  or_karo <- overlap/mean(unlist(resampling)[attr(unlist(resampling), "names") == "sample_overlap"])
  
  return(c(n_snps_overlap = overlap, or, or_perm = or_permutat, p_val = p_value, p_val_perm = p_value_permutat, p_val_emp = p_value_emp, p_value_karo = p_value_karo, or_karo = or_karo)) 
}
