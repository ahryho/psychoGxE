
source("integrative/util.R", chdir = TRUE)

pkg.list <- c("tidyverse", "dplyr", "data.table")
biocmanager.pkg.list <- c("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

LoadPackages(pkg.list)
LoadPackagesBio(biocmanager.pkg.list)

# Set up params

dmps.fn      <- '/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/methylation/01_lmem_dnam/dnam_lmem_svs_pcs_rslt.txt'
dmps.anno.fn <- '/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/methylation/02_dmp/dex_cpgs_annotated.csv'
 
# Load data
dmps.df <- fread(dmps.fn, sep = "\t")

# Prepare annotation tbl

anno.epic     <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.epic.sub <- as.data.frame(anno.epic[anno.epic$Name %in% dmps.df$PROBE_ID, ]) 

dmps.anno.df <- anno.epic.sub

# Save

write.table(dmps.anno.df, file = dmps.anno.fn, sep = "\t", row.names = F)
