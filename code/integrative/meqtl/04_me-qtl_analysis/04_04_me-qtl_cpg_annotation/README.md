# Functional enrichment analysis
Regarding functional enrichment in relation to nearby genes, the GR-response meQTL CpG sites were first overlapped with genomic annotation from UCSC for the hg19 genome build using TxDb.Hsapiens.UCSC.hg19.knownGene and ChIPseeker Bioconductor R packages21. Regarding functional enrichment in relation to CpG islands, the GR-induced meQTL CpG sites were mapped to the genome location according to Illumina’s annotation using minfi R package. After that, these cis-meQTL CpGs were tested for enrichment against baseline cis-meQTL CpGs.

# Histone mark enrichment analysis
To identify whether GR-response meQTL CpGs are enriched in specific chromatin states, the 20-state annotation inferred by the integrative and discriminative epigenome annotation system (IDEAS) from 18 blood and T- & B-cell cell tissue group (n = 18 cell lines) was used.  The position-based overlap of the GR-response meCpGs and chromatin states compared the overlap observed with 1000 equal-sized sets of baseline meCpGs.

### Folders structure

- `01_prepare the data`: 

### Results

The results are stored on the MPIP computational cluster: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls`