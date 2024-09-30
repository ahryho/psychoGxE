# Phenome-wide association studies (PheWAS) analysis

To investigate whether there is a causal relationship between the cell-type specific GR-meSNPs and psychiatric traits, a two-sample MR approach was employed using the R package TwoSampleMR.

_Genetic instruments_

The effect and variance estimates of identified SNP alleles were used as exposure phenotypes. The eight functional meQTL SNP groups were created based on the intersection between baseline and GR-response cell-type specific meQTL CpG sites significant at FDR < 0.05, and those meQTL CpGs significantly enriched in 20-states from IDEAS. These SNPs were grouped according to their associated CpG site (or multiple CpGs if applicable) and, to retain statistically independent SNPs as genetic instruments, clumped separately for each of eight SNP groups using a 10,000kB window and $R^2$ = 0.001 according to a European reference population. As a result, one single SNP per probe (a total of 154 SNPs) was retained for final analyses due to the high LD between SNPs associated with the same CpG site. 

_Outcome Phenotype_

As outcome phenotypes, the PGC cross-disorder phenotype obtained from a GWAS with one of eight psychiatric disorders (anorexia nervosa, ADHD, ASD, BD, MDD, obsessive-compulsive disorder, SCZ, and Tourette syndrome) was used in addition to 5,650 phenotypes originating from the UK Biobank study60 and the NHGRI-EBI GWAS Catalog. These phenotypes underwent the three-step filtration procedure described previously. The final list included 4,439 phenotypes combined with the cross-disorder phenotype.

_Mendelian Randomization_

Variant data from exposure and outcome phenotypes were harmonised by trying to infer positive strand alleles and removing ambiguous and palindromic SNPs. MR analyses were then applied to all remaining exposure-outcome pairs using Wald ratio estimation as a method of choice for single SNP MR. The obtained MR association estimates for all exposure-outcome combinations with available SNP data were corrected for multiple comparison corrections using the Benjamini-Hochberg FDR method.

To better interpret, all phenotypes with at least one FDR-significant association were grouped into the following seven categories: Biomarker, body morphology, disease, immune system, neuro-behavioural, nutrition, and other. A $χ^2$ test was applied for each phenotype category to assess deviations from the expected proportions of eight functional SNP groups.
The SNP-CpG pairs with most PheWAS associations in each phenotype category, and the outcomes that displayed associations with the largest number of SNP-CpGs were extracted.

## Folder structure

- `00_functions`: functions to extract and plot to snp hits and their categories
- `01_outcome_prep_ieu`: outcome selection using IEU GWAS Catalog. Script prepared by Nils Kappelmann
- `02_exposure_prep`: exposure data preparation
- `03_data_harmonisation`: data harmonisation
- `04_mr_analysis`: MR analysis
- `05_mr_result_processing`: process MR results obtained at the previous step
- `06a_mr_result_visualisation_8_gr`: MR results visualisation for 8 cell-type specific groups
- `06b_mr_result_visualisation_2_gr`: MR results visualisation for 2 cell-type specific groups

## Data

The exposure data, the IEU outcome list and trait annotation tables are store on the MPIP computational cluster.

## Results

The results are stored on the MPIP computational cluster.

## References

1. Zhang, Y., An, L., Yue, F. & Hardison, R. C. Jointly characterizing epigenetic dynamics across multiple human cell types. Nucleic Acids Res 44, 6721–6731 (2016).
2. Zhang, Y. & Hardison, R. C. Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation. Nucleic Acids Res 45, 9823–9836 (2017).
3. Lee, P. H. et al. Genomic Relationships, Novel Loci, and Pleiotropic Mechanisms across Eight Psychiatric Disorders. Cell 179, (2019).
60.	Bycroft, C. et al. The UK Biobank resource with deep phenotyping and genomic data. Nature 562, (2018).
61.	Buniello, A. et al. The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. Nucleic Acids Res 47, (2019).
62.	Penner-Goeke, S. et al. Assessment of glucocorticoid-induced enhancer activity of eSNP regions using STARR-seq reveals novel molecular mechanisms in psychiatric disorders. doi:10.1101/2022.05.18.22275090.
63.	Palmer, T. M. et al. Instrumental variable estimation of causal risk ratios and causal odds ratios in mendelian randomization analyses. Am J Epidemiol 173, (2011).

