library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation/EnrichmentWithPermutation_FUN.R")
source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation/RunEnrichmentChromHMM.R")

out.dir.pre             <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/"
# out.dir.pre           <- "~/bio/code/mpip/dex-stim-human-array/"
gwas.cluster.out.dir    <- paste0(out.dir.pre, "data/public_data/GWAS/")
gwas.enrichment.res.dir <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/with_missingness/enrichment/snps/")

meqtl.veh.snp.gr   <- readRDS(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtl_veh_snps_with_maf_gr.rds"))
meqtl.delta.snp.gr <- readRDS(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/meqtls/region_wise_independent_snps/meqtl_delta_snps_with_maf_gr.rds"))
background.all.gr  <- readRDS(paste0(out.dir.pre, "data/snps/imputed_qc/from_janine/qc/dex_geno_snps_with_maf_bins_gr.rds"))

# ChromHMM

chromhmm.blood.states <- readRDS(paste0(out.dir.pre, "data/annotation/chromHMM/chromHMM_blood_states.Rds"))
chromhmm.brain.states <- readRDS(paste0(out.dir.pre, "data/annotation/chromHMM/chromHMM_brain_states.Rds"))

# chromHMM Blood

## Delta vs Baseline

out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_blood_enrichment_delta_vs_veh.csv")

nperm  <- 10

chromhmm.enrich.perm.rslt <- RunEnrichmentChromHMM(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   chromhmm.states = chromhmm.blood.states,
                                                   out.fn = out.fn,
                                                   fun = EnrichmentWithPermutation,
                                                   nperm = nperm)

## Delta vs All

out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_blood_enrichment_delta_vs_all.csv")

nperm  <- 1000

chromhmm.enrich.perm.rslt <- RunEnrichmentChromHMM(own = meqtl.delta.snp.gr, 
                                                   background = background.all.gr, 
                                                   chromhmm.states = chromhmm.blood.states,
                                                   fun = EnrichmentWithPermutation,
                                                   out.fn = out.fn,
                                                   nperm = nperm)

## Baseline vs All

out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_blood_enrichment_veh_vs_all.csv")

nperm  <- 1000

chromhmm.enrich.perm.rslt <- RunEnrichmentChromHMM(own = meqtl.veh.snp.gr, 
                                                   background = background.all.gr, 
                                                   chromhmm.states = chromhmm.blood.states,
                                                   fun = EnrichmentWithPermutation,
                                                   out.fn = out.fn,
                                                   nperm = nperm)

# chromHMM Brain

## Delta vs Baseline

out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_brain_enrichment_delta_vs_veh.csv")

nperm  <- 1000

chromhmm.enrich.perm.rslt <- RunEnrichmentChromHMM(own = meqtl.delta.snp.gr, 
                                                   background = meqtl.veh.snp.gr, 
                                                   chromhmm.states = chromhmm.brain.states,
                                                   fun = EnrichmentWithPermutation,
                                                   out.fn = out.fn,
                                                   nperm = nperm)

## Delta vs All

out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_brain_enrichment_delta_vs_all.csv")

nperm  <- 1000

chromhmm.enrich.perm.rslt <- RunEnrichmentChromHMM(own = meqtl.delta.snp.gr, 
                                                   background = background.all.gr, 
                                                   chromhmm.states = chromhmm.brain.states,
                                                   fun = EnrichmentWithPermutation,
                                                   out.fn = out.fn,
                                                   nperm = nperm)

## Baseline vs All

out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_brain_enrichment_veh_vs_all.csv")

nperm  <- 1000

chromhmm.enrich.perm.rslt <- RunEnrichmentChromHMM(own = meqtl.veh.snp.gr, 
                                                   background = background.all.gr, 
                                                   chromhmm.states = chromhmm.brain.states,
                                                   fun = EnrichmentWithPermutation,
                                                   out.fn = out.fn,
                                                   nperm = nperm)

# Run ChromHMM run enrichment for each blood (brain) cell type separately 

# chromHMM Blood

tissue <- "blood"

if (!file.exists(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue)))  
  dir.create(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue))

## Delta vs Baseline

cell.types.lst <- elementMetadata(chromhmm.blood.states)[, "code"] %>% unique() %>% sort()

nperm  <- 1000

lapply(cell.types.lst, function(x){
  
  out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue, "/meqtl_snps_chromHMM_", tissue, "_enrichment_delta_vs_veh_", x, ".csv")
  
  RunEnrichmentChromHMM( 
       own = meqtl.delta.snp.gr, 
       background = meqtl.veh.snp.gr,
       chromhmm.states = chromhmm.blood.states[ elementMetadata(chromhmm.blood.states)[, "code"] == x, ],
       fun = EnrichmentWithPermutation,
       out.fn = out.fn,
       nperm = nperm)
  })

## Baseline vs all

tissue <- "blood"

out.dir.pre  <- "~/bio/code/mpip/dex-stim-human-array/"

if (!file.exists(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue, "_veh_vs_all")))  
  dir.create(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue, "_veh_vs_all"))

cell.types.lst <- elementMetadata(chromhmm.blood.states)[, "code"] %>% unique() %>% sort()

nperm  <- 1

lapply(cell.types.lst, function(x){
  
  out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue, "_veh_vs_all/meqtl_snps_chromHMM_", tissue, "_enrichment_veh_vs_all_", x, ".csv")
  
  RunEnrichmentChromHMM( 
    own = meqtl.veh.snp.gr, 
    background = background.all.gr,
    chromhmm.states = chromhmm.blood.states[ elementMetadata(chromhmm.blood.states)[, "code"] == x, ],
    fun = EnrichmentWithPermutation,
    out.fn = out.fn,
    nperm = nperm)
})

# chromHMM Brain

tissue <- "brain"

if (!file.exists(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue)))  
  dir.create(paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue))

## Delta vs Baseline

cell.types.lst <- elementMetadata(chromhmm.brain.states)[, "code"] %>% unique() %>% sort()

nperm  <- 1000

lapply(cell.types.lst[8:10], function(x){
  
  out.fn <- paste0(out.dir.pre, "output/data/integrative/matrixEQTL/05_me-qtl_enrichment/region_wise_independent_snps/meqtl_snps_chromHMM_", tissue, "/meqtl_snps_chromHMM_", tissue, "_enrichment_delta_vs_veh_", x, ".csv")
  
  RunEnrichmentChromHMM( 
    own = meqtl.delta.snp.gr, 
    background = meqtl.veh.snp.gr,
    chromhmm.states = chromhmm.blood.states[ elementMetadata(chromhmm.brain.states)[, "code"] == x, ],
    fun = EnrichmentWithPermutation,
    out.fn = out.fn,
    nperm = nperm)
})
