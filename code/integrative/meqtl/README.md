# Methylation quantitative trait loci analysis
The methylation QTL (meQTL) analysis was restricted to those SNP-CpG pairs that map within a region of 1Mb upstream or downstream. For each CpG-SNP pair, with the SNP encoded by 0, 1 and 2 according to the frequency of the minor allele and covariates, the association between baseline methylation site, and genotype, was assumed to be linear. To obtain baseline cis-meQTLs, the model was constructed controlling for sex, age, BMI, depression status (MDD), first three SVs derived from DNA methylation data, first two PCs which describes the genetic variability of the study cohort, and smoking score (SS) as described before and using the R package MatrixEQTL:

$CpG_{baseline}$ $∼$ β_0$  $+$  $β_1 \times$ Sex $+$ $β_2 \times$ Age $+$ $β_3 \times$ BMI $+$ $β_4 \times$ MDD $+$ $β_{5-7} \times SV_{5-7}$ $+$ $β_{8,9} \times PC_{1,2} $+$ $β_{10} \times SS$ $+$ $β_{SNP} \times$ SNP $+$ $ϵ$

To calculate GR-response cis-meQTLs, a linear model of the change on DNA methylation was constructed between baseline and GR-response using residuals from the linear regression and adjusting for the same covariates as for baseline meQTLs. 

$CpG_{GR}$ $=$ $res_{post-dex}$ $-$ $res_{baseline}$

$res_{post-dex}$ $∼$ $CpG_{post-dex}$ $-$ $(β_0$  $+$  $β_1 \times$ Sex $+$ $β_2 \times$ Age $+$ $β_3 \times$ BMI $+$ $β_4 \times$ MDD $+$ $β_{5-7} \times SV_{5-7}$ $+$ $β_{8,9} \times PC_{1,2}$ $+$ $β_{10} \times SS$ $+$ $β_{SNP} \times$ SNP $+$ $ϵ)$

$res_{basline}$ $∼$ $CpG_{baseline}$ $-$ $(β_0$  $+$  $β_1 \times$ Sex $+$ $β_2 \times$ Age $+$ $β_3 \times$ BMI $+$ $β_4 \times$ MDD $+$ $β_{5-7} \times SV_{5-7}$ $+$ $β_{8,9} \times PC_{1,2}$ $+$ $β_{10} \times SS$ $+$ $β_{SNP} \times$ SNP $+$ $ϵ)$

Since cis-regions with an extensive LD structure increase the false positive meQTLs, the Benjamini-Hochberg FDR method was applied to correct the adjusted p-value significance by using only the most significant and independent SNPs per probe. The number of independent meSNPs per cis region was identified by LD clumping the SNPs using the clump command in PLINK. Each independent SNP formed a SNP bin aggregating all other SNPs into bins by independent SNP at $r_2$ $>$ 0.2 and distance < 1Mb, such that all SNPs within a given bin were correlated to the independent SNP, but to any other tag SNP.  The false-positive SNP-probe pairs were limited to less than 5%, applying the **FDR** as statistical significance.

### Folders structure

- `01_prepare the data`: 
  - `dnam_residuals`: calculation of DNAm residuals
  - `01_prepare_data`: script for preparing omic data for QTL analysis, i.e. input for MatrixEQTL 
  - `02_prepare_meqtl_parallel_and_opposite_fc_groups`: splitting the meqtl into groups based on allelic direction
  - `03_prepare_snps`: script for preparing SNP data for MatrixEQTL
- `02_get_distances_between_cpgs_ensg`: calculation distances between CpGs and ENSGs
- `03_run_qtl_analysis`: meQTLs calculation using MatrixEQTL package
- `04_me-qtl_analysis`: 
  - `04_01_get_clumped_meqtl_rw`: region-wise clumping of meQTLs
  - `04_02_most_likely_genotype`: analysis of meQTL associations using SNPs in which missing values got replaced with the most likely genotype. The folder contains the analysis of both primary and clumped associations, as well as functional enrichment, chromHMM, GWAS enrichemnt analyses of meQTL CpGs and meQTL SNPs. In addition, it also contains analysis of cell-type specific meQTL CpGs
  - `04_04_me-qtl_cpg_annotation`: functional genomic and histone mark annotation and enrichment analyses of meQTL CpGs. For more details, please refer to the [folder](https://github.com/ahryho/psychoGE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_cpg_annotation).
  - `04_04_me-qtl_snp_annotation`: functional genomic and histone mark annotation and enrichment analyses, GWAS enrichment analysis of meQTL SNPs. For more details, please refer to the [folder](https://github.com/ahryho/psychoGE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation).
  - `04_05_overlap_mesnp_with_public_data`: comparison of meQTL SNPs to GR-eQTL SNPs from Moore et al.
  - `04_06_with_missigness`: analysis of meQTL associations using SNPs with missing genotypes. The folder contains the analysis of both primary and clumped associations, as well as functional enrichment, chromHMM, GWAS enrichemnt analyses of meQTL CpGs and meQTL SNPs. In addition, it also contains analysis of cell-type specific meQTL CpGs
  - `04_07_comparison_between_meqtl_approaches`: comparison between meQTL obtained by using different sets of SNPs: (a) with missing values, (b) probablities, (c) most linkely genotype
  - `04_08_allelic_direction_groups`: analysis of meQTL groups. Groups were created based on the genotype allelic direction.

### Results

The results are stored on the MPIP computational cluster.

