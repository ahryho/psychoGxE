library(data.table)
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)

library(parallel)
library(foreach)
library(doParallel)

source("~/mpip/projects/dex-stim-human-array/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_cpg_annotation/EnrichmentWithPermutation_FUN.R")

out.dir.pre  <- "/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/"
setwd(out.dir.pre)

# Relation to island

meqtls.cpg.delta.coord.gr <- readRDS("annotation/cpgs/meqtl_delta_cpgs_gr.rds")
meqtls.cpg.veh.coord.gr   <- readRDS("annotation/cpgs/meqtl_veh_cpgs_gr.rds")

all.cpg.coord.gr          <- readRDS("annotation/cpgs/all_cpgs_gr.rds")
anno.epic.gr              <- readRDS("/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/annotation/anno_epic_gr.rds")

features.lst <- elementMetadata(meqtls.cpg.delta.coord.gr)[, "Relation_to_Island"] %>% unique() %>% sort()

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                       state <- features.lst[i]                                     
                                       public <- anno.epic.gr[(elementMetadata(anno.epic.gr)[, "Relation_to_Island"]) == state, ] 
                                       EnrichmentWithPermutationGeneLocWithoutMAF(own = meqtls.cpg.delta.coord.gr,
                                                                           background = meqtls.cpg.veh.coord.gr,
                                                                           public = public,
                                                                           nperm = nperm)
                                     }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Feature = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

fwrite(gen.loc.enrich.perm.rslt, 
       file = "enrichment/cpgs/meqtl_cpgs_relation_to_island_enrichment_perm.csv", 
       row.names = F, quote = F)

# Relation to island

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                      state <- features.lst[i]                                     
                                      public <- anno.epic.gr[(elementMetadata(anno.epic.gr)[, "Relation_to_Island"]) == state, ] 
                                      EnrichmentWithPermutationGeneLocWithoutMAF(own = meqtls.cpg.delta.coord.gr,
                                                                                 background = all.cpg.coord.gr,
                                                                                 public = public,
                                                                                 nperm = nperm)
                                    }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Feature = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

fwrite(gen.loc.enrich.perm.rslt, 
       file = "enrichment/cpgs/meqtl_cpgs_relation_to_island_enrichment_perm_delta_vs_all.csv", 
       row.names = F, quote = F)

# Relation to island

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                      state <- features.lst[i]                                     
                                      public <- anno.epic.gr[(elementMetadata(anno.epic.gr)[, "Relation_to_Island"]) == state, ] 
                                      EnrichmentWithPermutationGeneLocWithoutMAF(own = meqtls.cpg.veh.coord.gr,
                                                                                 background = all.cpg.coord.gr,
                                                                                 public = public,
                                                                                 nperm = nperm)
                                    }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Feature = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

fwrite(gen.loc.enrich.perm.rslt, 
       file = "enrichment/cpgs/meqtl_cpgs_relation_to_island_enrichment_perm_veh_vs_all.csv", 
       row.names = F, quote = F)


###
# ChIPSeeker 
###

delta.meqtl.cpg.anno.gr <- readRDS("annotation/cpgs/meqtls_cpg_annotated_withChIPseeker_delta_gr.rds")
veh.meqtl.cpg.anno.gr   <- readRDS("annotation/cpgs/meqtls_cpg_annotated_withChIPseeker_veh_gr.rds")
all.cpg.anno.gr         <- readRDS("annotation/cpgs/cpg_annotated_withChIPseeker_all_gr.rds")

features.lst <- elementMetadata(delta.meqtl.cpg.anno.gr)[, "annotation"] %>% unique() %>% sort()

# Delta vs Baseline

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                      feature <- features.lst[i]                                     
                                      EnrichmentWithPermutationWithoutMAFnoPubData(own = delta.meqtl.cpg.anno.gr,
                                                                                   background = veh.meqtl.cpg.anno.gr,
                                                                                   feature = feature,
                                                                                   nperm = nperm)
                                    }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Feature = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

fwrite(gen.loc.enrich.perm.rslt, 
       file = "enrichment/cpgs/meqtl_cpgs_chipseeker_enrichment_perm_delta_vs_veh.csv", 
       row.names = F, quote = F)

# Delta vs All

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                      feature <- features.lst[i]                                     
                                      EnrichmentWithPermutationWithoutMAFnoPubData(own = delta.meqtl.cpg.anno.gr,
                                                                                   background = all.cpg.anno.gr,
                                                                                   feature = feature,
                                                                                   nperm = nperm)
                                    }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Feature = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

fwrite(gen.loc.enrich.perm.rslt, 
       file = "enrichment/cpgs/meqtl_cpgs_chipseeker_enrichment_perm_delta_vs_all.csv", 
       row.names = F, quote = F)

# Baseline vs All

no.cores <- detectCores() - 2
cl <- makeCluster(no.cores)
registerDoParallel(cl)

nperm   <- 1000

gen.loc.enrich.perm.rslt <- foreach(i =  seq_along(features.lst), 
                                    .combine = rbind, 
                                    .packages =  c("GenomicRanges", "dplyr")) %dopar% {
                                      feature <- features.lst[i]                                     
                                      EnrichmentWithPermutationWithoutMAFnoPubData(own = veh.meqtl.cpg.anno.gr,
                                                                                   background = all.cpg.anno.gr,
                                                                                   feature = feature,
                                                                                   nperm = nperm)
                                    }

stopImplicitCluster()

gen.loc.enrich.perm.rslt <- cbind(gen.loc.enrich.perm.rslt %>% data.frame(row.names = NULL), 
                                  Feature = features.lst)

gen.loc.enrich.perm.rslt[["n_perm"]] <- nperm

gen.loc.enrich.perm.rslt

fwrite(gen.loc.enrich.perm.rslt, 
       file = "enrichment/cpgs/meqtl_cpgs_chipseeker_enrichment_perm_veh_vs_all.csv", 
       row.names = F, quote = F)
