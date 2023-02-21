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

# Run ChromHMM run enrichment for each blood (brain) cell type separately 

# chromHMM Blood

tissue <- "blood"

if (!file.exists(paste0("enrichment/cpgs/meqtl_cpgs_chromHMM_", tissue)))
  dir.create(paste0("enrichment/cpgs/meqtl_cpgs_chromHMM_", tissue))

## Load ChromHMM states

chromhmm.blood.states <- readRDS("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/annotation/chromHMM/chromHMM_blood_states.Rds")

## Delta vs Baseline

cell.types.lst <- elementMetadata(chromhmm.blood.states)[, "code"] %>% unique() %>% sort()

nperm   <- 10000 

lapply(cell.types.lst, function(x){
  
  out.fn <- paste0("enrichment/cpgs/meqtl_cpgs_chromHMM_", tissue, "/meqtl_cpgs_chromHMM_", tissue, "_enrichment_delta_vs_veh_", x, ".csv")
  
  RunEnrichmentChromHMM( 
    own = meqtl.delta.gr, 
    background = meqtl.veh.gr,
    chromhmm.states = chromhmm.blood.states[ elementMetadata(chromhmm.blood.states)[, "code"] == x, ],
    fun = EnrichmentWithPermutationGeneLocWithoutMAF,
    out.fn = out.fn,
    nperm = nperm)
})


# chromHMM Brain

tissue <- "brain"

if (!file.exists(paste0("enrichment/cpgs/meqtl_cpgs_chromHMM_", tissue)))
  dir.create(paste0("enrichment/cpgs/meqtl_cpgs_chromHMM_", tissue))

## Load ChromHMM states

chromhmm.brain.states <- readRDS("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/annotation/chromHMM/chromHMM_brain_states.Rds")

## Delta vs Baseline

cell.types.lst <- elementMetadata(chromhmm.brain.states)[, "code"] %>% unique() %>% sort()

nperm   <- 10000 

lapply(cell.types.lst, function(x){
  
  out.fn <- paste0("enrichment/cpgs/meqtl_cpgs_chromHMM_", tissue, "/meqtl_cpgs_chromHMM_", tissue, "_enrichment_delta_vs_veh_", x, ".csv")
  
  RunEnrichmentChromHMM( 
    own = meqtl.delta.gr, 
    background = meqtl.veh.gr,
    chromhmm.states = chromhmm.brain.states[ elementMetadata(chromhmm.brain.states)[, "code"] == x, ],
    fun = EnrichmentWithPermutationGeneLocWithoutMAF,
    out.fn = out.fn,
    nperm = nperm)
})
