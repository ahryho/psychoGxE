# GWAS and histone mark enrichemnts of meQTL SNPs

## Histone mark enrichment analysis

To identify whether GR-response meQTL SNPs are enriched in specific chromatin states, the 20-state annotation inferred by the integrative and discriminative epigenome annotation system (IDEAS) from 18 blood and T- & B-cell cell tissue group (n = 18 cell lines) was used.  The position-based overlap of the GR-response meQTL SNPs and chromatin states compared the overlap observed with 1000 equal-sized sets of baseline meQTL SNPs. The enrichment analysis was carried out adjusting for MAF.

## GWAS enrichment analysis

Concerning enrichment for GWAS hits, the independent GR-response cis-meQTL SNPs were matched to GWAS variants based on chromosome and position (hg19). The full summary statistics were derived from the respective publication to check for enrichment for nominal significant GWAS hits with regard to the independent baseline cis-meQTL SNPs. The enrichment meQTL SNPs and GWAS risk-SNPs were tested compared to 1000 MAF-matched independent meSNP sets. The Fisher-test empirical p-value $\leq$ 0.05 was used to define the significant differences. A total of 1,000 permutations were carried out per GWAS.

## Folders structure

- `01_prepare the data`: annotate meQTL SNPs, prepare data for the enrichemnt analysis
- `02a_generate_gr_functions`: function to generate GRanges objects for enrichment 
- `02b_generate_gr`: generate GRanges of meQTL SNPs and GWASes
- `03a_run_GWAS_enrichment`: GWAS enrichment of meQTL SNPs
- `03b_run_chromHMM_enrichment`: histone mark enrichment of meQTL SNPs
- `EnrichmentWithPermutation_FUN`: function for the enrichment with permutation
- `RunEnrichmentChromHMM`: fucntion for the histone mark enrichment
- `PlotChromHMMEnrichment`: function for visualisation of histone mark enrichment results
- `SavePPTXChromHMMEnrichmentPlot`: Function to save enrichment plots to .pptx file

## Results

The results are stored on the MPIP computational cluster.

## References

1. Zhang, Y., An, L., Yue, F. & Hardison, R. C. Jointly characterizing epigenetic dynamics across multiple human cell types. Nucleic Acids Res 44, 6721–6731 (2016).
2. Zhang, Y. & Hardison, R. C. Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation. Nucleic Acids Res 45, 9823–9836 (2017).
3. Grove, J. et al. Mette Nyegaard 1,2,3 , Terje Naerland 14,30 , Duncan S. Palmer 5,6 , Aarno Palotie 5,6,22,31 , Carsten Bøcker Pedersen 1,12,13 , Marianne Giørtz Pedersen 1,12,13 , Timothy dPoterba 5. Autism Spectrum Disorder Working Group of the Psychiatric Genomics Consortium 23, 22.
4. Mullins, N. et al. Genome-wide association study of over 40,000 bipolar disorder cases provides new insights into the underlying biology. Nat Genet 53, 817 (2021).
5. Giannakopoulou, O. et al. The Genetic Architecture of Depression in Individuals of East Asian Ancestry: A Genome-Wide Association Study. JAMA Psychiatry 78, 1258–1269 (2021).
6. Trubetskoy, V. et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. Nature 604, 502–508 (2022).
7. Coleman, J. R. I. et al. The Genetics of the Mood Disorder Spectrum: Genome-wide Association Analyses of More Than 185,000 Cases and 439,000 Controls. Biol Psychiatry 88, 169–184 (2020).

