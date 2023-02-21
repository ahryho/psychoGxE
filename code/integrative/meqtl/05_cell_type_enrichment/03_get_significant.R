require(arrow)
require(data.table)

GetFDRdf <- function(df, treatment, out.fn){
  fdr.df <- apply(df[, -1], 2, p.adjust, method = "fdr")
  fdr.df <- data.frame(cbind(CpG_ID = df$CpG_ID, fdr.df))

  reshaped.fdr.df <- fdr.df %>% 
    reshape2::melt(measure.vars = colnames(fdr.df)[2:13]) %>% setDT()
  colnames(reshaped.fdr.df) <- c("CpG_ID", "Type", "fdr")
  
  reshaped.fdr.df <- reshaped.fdr.df[fdr <= 0.05]
  reshaped.fdr.df[["Treatment"]] <- treatment
  
  fwrite(reshaped.fdr.df,
         out.fn,
         quote = F, row.names = F, sep = ";")
}

# Subset baseline and dex

bcc.rslt.dir <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/cell_type_enrichment/"
setwd(bcc.rslt.dir)

for (treatment in c("dex", "veh")){
  bcc.df <- read_delim_arrow(file = paste0("dnam_cell_type_enrichment_pvals_", treatment, ".csv"), 
                             delim = ";",
                             col_select = 1:13)

  out.fn <- paste0("dnam_cell_type_enrichment_", treatment, "_fdr_005.csv")

  GetFDRdf(bcc.df, treatment, out.fn)
}

# Subset meqtl CpGs from baseline and dex

dex.bcc.fdr.df <- fread(paste0(bcc.rslt.dir, "dnam_cell_type_enrichment_dex_fdr_005.csv"))
veh.bcc.fdr.df <- fread(paste0(bcc.rslt.dir, "dnam_cell_type_enrichment_veh_fdr_005.csv"))

## load delta  meqtls
out.dir.pre      <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/"
meqtl.delta.df   <- fread(paste0(out.dir.pre, "clumped/me-qtl_cis_result_snps_with_na_delta_beta_fdr_005.csv"), col.names = c("SNP", "CpG_ID", "beta", "t-stat", "p-value", "fdr"))
delta.meqtl.cpgs <- meqtl.delta.df$CpG_ID %>% unique() 

## save the subset 
out.fn            <- paste0(out.dir.pre, "enrichment/meqtl_veh_cell_type_enrichment_fdr_005.csv")
fwrite(veh.bcc.fdr.df[CpG_ID %in% delta.meqtl.cpgs, ], out.fn,
       quote = F, row.names = F, sep = ";")

out.fn            <- paste0(out.dir.pre, "enrichment/meqtl_dex_cell_type_enrichment_fdr_005.csv")
fwrite(dex.bcc.fdr.df[CpG_ID %in% delta.meqtl.cpgs, ], out.fn,
       quote = F, row.names = F, sep = ";")
