# Methylation quantitative trait loci analysis
The methylation QTL (meQTL) analysis was restricted to those SNP-CpG pairs that map within a region of 1Mb upstream or downstream. For each CpG-SNP pair, with the SNP encoded by 0, 1 and 2 according to the frequency of the minor allele and covariates, the association between baseline methylation site, and genotype, was assumed to be linear. To obtain baseline cis-meQTLs, the model was constructed controlling for sex, age, BMI, depression status (MDD), first three SVs derived from DNA methylation data, first two PCs which describes the genetic variability of the study cohort, and smoking score (SS) as described before and using the R package MatrixEQTL:

$CpG_{baseline} ∼ β_0$  $+$  $β_1 \times$ Sex $+$ $β_2 \times$ Age $+$ $β_3 \times$ BMI $+$ $β_4 \times$ MDD $+$ $β_(5-7) \times SV_{5-7} $+$ $β_{8,9} \times PC_{1,2} $+$ $β_{10} \times SS$ $+$ $β_{SNP} \times$ SNP $+$ $ϵ$

To calculate GR-response cis-meQTLs, a linear model of the change on DNA methylation was constructed between baseline and GR-response using residuals from the linear regression and adjusting for the same covariates as for baseline meQTLs. 

$CpG_{GR} ∼ $β_0$ $+$ $β_{SNP} \times$ SNP $+$ $ϵ$

CpG_GR= res_(post-dex) " - " res_baseline

res_(post-dex)∼〖CpG_(post-dex)+ β〗_0 " + " β_1×"Sex"+β_2×"Age"+β_3×"BMI"+β_4×"MDD"+β_(5-7)×〖"SV" 〗_(5-7)+β_8,9×〖"PC" 〗_8,9+ β_10×SS" + " β_SNP×"SNP"+ϵ

res_baseline∼〖CpG_baseline+ β〗_0 " + " β_1×"Sex"+β_2×"Age"+β_3×"BMI"+β_4×"MDD"+β_(5-7)×〖"SV" 〗_(5-7)+β_8,9×〖"PC" 〗_1,2+ β_10×SS" + " β_SNP×"SNP"+ϵ

Since cis-regions with an extensive LD structure increase the false positive meQTLs, the Benjamini-Hochberg FDR method was applied to correct the adjusted p-value significance by using only the most significant and independent SNPs per probe. The number of independent meSNPs per cis region was identified by LD clumping the SNPs using the clump command in PLINK. Each independent SNP formed a SNP bin aggregating all other SNPs into bins by independent SNP at $r_2 > 0.2$ and distance < 1Mb, such that all SNPs within a given bin were correlated to the independent SNP, but to any other tag SNP.  The false-positive SNP-probe pairs were limited to less than 5%, applying the **FDR** $< 5%$ as statistical significance.


### Folders structure

- `00_prepare the data`: the script for phenotype data to identify differentially methylated regions (DMRs)
- `01_tca`: Tensor Composition Analysis usig the R package TCA
- `02_lm_vif`: script for checking for multicollinarity between blood cell types using VIF
- `02a_lm`: cell type-specificity on DMAm using lm excluding the cell type with the highest VIF
- `02b_lm`: cell type-specificity on DMAm using lm including all cell types
- `03_get_significant`: get statistically signififcant at FDR < 0.05 associations, and extract assocaition with meQTL CpGs
- `04_epistress_score`: epistress score calculation for Sarah Merrill from the University of British Columbia 

### Results

The results are stored on the MPIP computational cluster: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/cell_type_enrichment/`