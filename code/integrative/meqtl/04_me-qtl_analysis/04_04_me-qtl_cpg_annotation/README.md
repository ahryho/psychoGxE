# Functional genomic and histone mark enrichemnts of meQTL CpGs

## Functional enrichment analysis

Regarding functional enrichment in relation to nearby genes, the GR-response meQTL CpG sites were first overlapped with genomic annotation from UCSC for the hg19 genome build using TxDb.Hsapiens.UCSC.hg19.knownGene and ChIPseeker Bioconductor R packages21. Regarding functional enrichment in relation to CpG islands, the GR-induced meQTL CpG sites were mapped to the genome location according to Illumina’s annotation using minfi R package. After that, these cis-meQTL CpGs were tested for enrichment against baseline cis-meQTL CpGs.

## Histone mark enrichment analysis

To identify whether GR-response meQTL CpGs are enriched in specific chromatin states, the 20-state annotation inferred by the integrative and discriminative epigenome annotation system $(IDEAS)^{1,2}$ from 18 blood and T- & B-cell cell tissue group (n = 18 cell lines) was used.  The position-based overlap of the GR-response meCpGs and chromatin states compared the overlap observed with 1000 equal-sized sets of baseline meCpGs.

## Folders structure

- `01_prepare the data`: annotate meQTL CpGs, prepare data for the enrichemnt analysis
- `02a_run_enrichment`: funcional genomic enrichment of meQTL CpGs
- `02b_run_enrichment_chromHMM`: histone mark enrichment of meQTL CpGs using Roadmap chromHMM
- `02c_run_enrichment_IDEAS`: histone mark enrichment of meQTL CpGs using IDEAS
- `04_me-qtl_cpg_annotation_rw_report`: visualisation of the results
- `05_summary_all_sign_enrich_terms`: summary of the results
- `EnrichmentWithPermutation_FUN`: functions for enrichment

## Results

The results are stored on the MPIP computational cluster: 

- annotation: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/annotation/cpgs`
- enrichment: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with missigness/enrichment/cpgs`

Visualisation of the results is [here](https://github.com/ahryho/psychoGE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_06_snps_with_missingness).

## References

1. Zhang, Y., An, L., Yue, F. & Hardison, R. C. Jointly characterizing epigenetic dynamics across multiple human cell types. Nucleic Acids Res 44, 6721–6731 (2016).
2. Zhang, Y. & Hardison, R. C. Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation. Nucleic Acids Res 45, 9823–9836 (2017).
