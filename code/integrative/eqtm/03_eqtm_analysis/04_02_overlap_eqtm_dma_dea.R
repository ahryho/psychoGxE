library(data.table)
library(arrow)
library(dplyr)
library(ggplot2)
library(eulerr)

# Set up global params

fc <- 0.02
fdr <- 0.01

# Load data
out.dir.pre <- "~/bio/code/mpip/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms/"

## map tbl
map.ilmn.ensg.gene.df <- read_delim_arrow("~/bio/code/mpip/dex-stim-human-array/data/mapping/mapping_ilmn_ensg_gene.csv", delim = ";") %>% setDT()
map.cpg.gene.ensg.df <- read_delim_arrow("~/bio/code/mpip/dex-stim-human-array/data/mapping/mapping_cpg_gene_ensg_full.csv", delim = ";") %>% setDT()

## eQTMs
col.names <- c("CpG_ID", "ENSG_ID", "beta", "t-stat", "p-value", "FDR")
sign.eqtm.df <- read_delim_arrow(paste0(out.dir.pre, "eqtm_cis_result_delta_no_residuals_or_cov_beta.csv"), delim = ";") %>% setDT()

sign.eqtm.df <- sign.eqtm.df[FDR < 0.05]

# sign.eqtm.df <- left_join(sign.eqtm.df, map.df, by = c("ENSG_ID" = "Ensemble_ID")) 
# sign.eqtm.df <- sign.eqtm.df %>% data.table()

## Significnat differential transcripts 

sign.gex.df <- read_delim_arrow("~/bio/code/mpip/dex-stim-human-array/output/data/gene_expression/01_lmem_dea/lmem_dex_svs_1_5.csv", delim = "\t") %>% setDT()

sign.gex.df <- sign.gex.df[pFDR < fdr][abs(FC) > fc]

sign.gex.df[FC > fc, reg := "UP"]
sign.gex.df[FC < fc, reg := "DOWN"]

### map ILMN_ID to ENSG_ID
gex.ensg.ids <- map.ilmn.ensg.gene.df[Illumina_ID %in% sign.gex.df$Probe_Id, Ensemble_ID] %>% unique() # 9'165

## Significnat dDMPs (FC > 0.02, FDR < 0.01)

sign.dmp.df <- fread("~/bio/code/mpip/dex-stim-human-array/output/data/methylation/01_lmem_dnam/dnam_lmem_svs_pcs_rslt.txt")[FDR < fdr][abs(FC) > fc]

sign.dmp.df[FC > fc, reg := "hyper"]
sign.dmp.df[FC < fc, reg := "hypo"]

### map CpG_ID to ENSG_ID
meth.ensg.ids <- map.cpg.gene.ensg.df[CpG_ID %in% sign.dmp.df$PROBE_ID, Ensemble_ID ] %>% unique() # 4'179

# Overlap

fit <- euler(list(DMP = meth.ensg.ids, DEA = gex.ensg.ids))
round(fit$original.value["DMP&DEG"] / length(gex.ensg.ids) * 100, 1)
round(fit$original.value["DMP&DEG"] / length(meth.ensg.ids) * 100, 1)

plot(fit,
     quantities = list(type = c("counts", "percent")),
     fill = c("red", "white"))

length(sign.gex.df$Probe_Id); length(gex.ensg.ids)
length(sign.dmp.df$PROBE_ID); length(meth.ensg.ids)


# Overlap 2
assign(gene.type, "Gene_ID")

gex.gene.up.ids <- map.ilmn.ensg.gene.df[Illumina_ID %in% sign.gex.df[reg == "UP", Probe_Id], get(gene.type)] %>% unique()
length(sign.gex.df[reg == "UP", Probe_Id]); length(gex.gene.up.ids); # 4466; 4486

gex.gene.down.ids <- map.ilmn.ensg.gene.df[Illumina_ID %in% sign.gex.df[reg == "DOWN", Probe_Id], get(gene.type) ] %>% unique()
length(sign.gex.df[reg == "DOWN", Probe_Id]); length(gex.gene.down.ids); # 4526; 4794

meth.gene.hyper.ids <- map.cpg.gene.ensg.df[CpG_ID %in% sign.dmp.df[reg == "hyper", PROBE_ID], get(gene.type) ] %>% unique()
length(sign.dmp.df[reg == "hyper", PROBE_ID]); length(meth.gene.hyper.ids); # 2072; 1165

meth.gene.hypo.ids <- map.cpg.gene.ensg.df[CpG_ID %in% sign.dmp.df[reg == "hypo", PROBE_ID], get(gene.type) ] %>% unique()
length(sign.dmp.df[reg == "hypo", PROBE_ID]); length(meth.gene.hypo.ids); # 7818; 3435

fit2 <- euler(list(dea_up = gex.gene.up.ids,
                   dea_down = gex.gene.down.ids,
                   dmp_hypo = meth.gene.hypo.ids, 
                   dmp_hyper = meth.gene.hyper.ids
))

plot(fit2,
     quantities = list(type = c("counts", "percent")),
     fills = palette.colors(palette = "Pastel 2")[6:9],
     labels = c("UP", "DOWN", "HYPO", "HYPER")
)

## Map.df for significant gex

map.sign.df <- map.df[Illumina_ID %in% sign.gex.df$Probe_Id, .(Illumina_ID, Gene_ID)]
map.sign.df <- map.sign.df[!duplicated(map.sign.df)]

# Prepare data for overlap
up.dea  <- sign.gex.df[FC > fc] # 4'466
down.dea <- sign.gex.df[FC < fc] # 4'526

hyper.dmp <- sign.dmp.df[FC > fc] # 2'072
hypo.dmp <- sign.dmp.df[FC < fc] # 7'818

sign.eqtm.df$CpG_ID %>% unique() %>% length() # 34'522
sign.eqtm.df$ENSG_ID %>% unique() %>% length() # 7'291

# Check how many hypo and hyper DMPs in eQTMs
table(hypo.dmp$PROBE_ID %in% sign.eqtm.df$CpG_ID) # 5441
table(hyper.dmp$PROBE_ID %in% sign.eqtm.df$CpG_ID) # 929

eqtms.hypo.dmps <- sign.eqtm.df[CpG_ID %in% hypo.dmp$PROBE_ID] # 25072
eqtms.hyper.dmps <- sign.eqtm.df[CpG_ID %in% hyper.dmp$PROBE_ID] # 4289

# Check how many up and down regulated gex in eQTMs with significant hypomethylated CpGs
eqtms.hypo.dmps$ENSG_ID %>% unique() %>% length() # 3'682
table(up.dea$Probe_Id %in% eqtms.hypo.dmps$Illumina_ID) # 2'302
table(down.dea$Probe_Id %in% eqtms.hypo.dmps$Illumina_ID) # 2'320