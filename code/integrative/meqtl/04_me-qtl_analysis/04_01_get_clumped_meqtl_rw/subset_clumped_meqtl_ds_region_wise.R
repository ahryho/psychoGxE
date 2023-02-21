library(data.table)
library(dplyr)

#' Function to subset clumped meqtls from meqtl dataste which contains statistics
#'
#' @param meqtl.fn 
#' @param clumped.meqtl.lst.fn 
#' @param out.fn 
#'
#' @return
#' @export
#'
#' @examples
SubsetClumpedQTLs <- function(meqtl.fn, clumped.meqtl.lst.fn, out.fn){
  meqtl.df          <- fread(meqtl.fn, col.names = c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr"))
  clumped.meqtl.lst <- fread(clumped.meqtl.lst.fn, header = F, col.names = c("CpG_ID", "SNP"))
  meqtl.clumped.df  <- left_join(clumped.meqtl.lst, meqtl.df)[,  c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr")]
  
  fwrite(meqtl.clumped.df,
         out.fn,
         quote = F, row.names = F, sep = ";")
}

# Subset clumped meqtls for each timepoint 

## Genotype probabilities

out.dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/snp_probabilities/"
setwd(out.dir.pre)

for (treatment in c("delta", "veh")){
  meqtl.fn             <- paste0("me-qtl_cis_result_snps_", treatment, "_probabilities_beta_fdr_005.csv")
  clumped.meqtl.lst.fn <- paste0("clumped/me-qtl_cis_", treatment, "_snp_probabilities_beta_clumped_cpg_snp_associations.txt")
  out.fn               <- paste0("clumped/me-qtl_cis_result_", treatment, "_snp_probabilities_beta_fdr_005.csv")
  
  SubsetClumpedQTLs(meqtl.fn, clumped.meqtl.lst.fn, out.fn)
}

## SNPs with NAs

out.dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/"
setwd(out.dir.pre)

for (treatment in c("delta", "veh")){
  meqtl.fn             <- paste0("me-qtl_cis_result_snps_with_na_",  treatment, "_beta_fdr_005.csv")
  clumped.meqtl.lst.fn <- paste0("clumped/me-qtl_cis_", treatment, "_snps_with_na_beta_clumped_cpg_snp_associations.txt")
  out.fn               <- paste0("clumped/me-qtl_cis_result_snps_with_na_", treatment, "_beta_fdr_005.csv")
  
  SubsetClumpedQTLs(meqtl.fn, clumped.meqtl.lst.fn, out.fn)
}

