library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)
library(plyranges)
library(BRGenomics)

library(parallel)
library(foreach)
library(doParallel)

out.dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/"
setwd(out.dir.pre)

source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_cpg_annotation/EnrichmentWithPermutation_FUN.R")
source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation/RunEnrichmentChromHMM.R")

meqtl.delta.gr <- readRDS("annotation/cpgs/meqtl_delta_cpgs_gr.rds")
meqtl.veh.gr   <- readRDS("annotation/cpgs/meqtl_veh_cpgs_gr.rds")

# Run ENCODE IDEAS run enrichment for each blood cell type separately 

# IDEAS Blood

tissue <- "blood"

if (!file.exists(paste0("enrichment/cpgs/meqtl_cpgs_IDEAS_", tissue)))
  dir.create(paste0("enrichment/cpgs/meqtl_cpgs_IDEAS_", tissue))

## Load IDEAS states

ideas.blood.states <- readRDS("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/public_data/ENCODE/IDEAS/all.blood.tissues.ideas.states.rds")
names(elementMetadata(ideas.blood.states))[1] <- "type"

## Delta vs Baseline

cell.types.lst <- elementMetadata(ideas.blood.states)[, "code"] %>% unique() %>% sort()

nperm   <- 10000

lapply(cell.types.lst[-c(1:8)], function(x){
  
  out.fn <- paste0("enrichment/cpgs/meqtl_cpgs_IDEAS_", tissue, "/meqtl_cpgs_IDEAS_", tissue, "_enrichment_delta_vs_veh_", x, ".csv")
  
  RunEnrichmentChromHMM( 
    own = meqtl.delta.gr, 
    background = meqtl.veh.gr,
    chromhmm.states = ideas.blood.states[ elementMetadata(ideas.blood.states)[, "code"] == x, ],
    fun = EnrichmentWithPermutationGeneLocWithoutMAF,
    out.fn = out.fn,
    nperm = nperm)
})

