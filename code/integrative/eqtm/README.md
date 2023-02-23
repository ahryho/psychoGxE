## Expression quantitative trait methylation analysis

The baseline/GR-response expression quantitative trait methylation (eQTM) analysis was conducted on 11,944 baseline/GR-response gene expressions and baseline/GR-response 740,357 CpG sites using linear regression. Each model was adjusted for covariates: sex, age, BMI, depression status (MDD), first three SVs derived from DNA methylation data, first two PCs from genotypic data and Smoking score (SS):  

GEX $∼$ $β_0 + β_{CpG} \times$ CpG $+$ $β_1 \times$ Sex $+$ $β_2 \times$ Age $+$ $β_3 \times$ BMI $+$ $β4_4 \times$ MDD $+$ $β_{5-7} \times SV_{5-7}$ $+$ $β_{8,9} \times PC_{1,2}$ $+$ $β_{10} \times$ SS $+$ $ϵ$ 

The eQTM analysis aims to investigate whether the DNA methylation influences gene expression and the integration of the results with meQTL analysis. Therefore, the nominal p-value of 0.05 was used to identify the number of significant baseline and GR-response eQTM pairs.

### Folders structure

- `01_prepare the data`: the script for correcting gene-expression data for covariates
- `02_run_eqtm`: the scripts for DMPs identifications. 6 models were tested:
  
    - `run_eqtm`: eQTMs calculation using MatrixEQTL package
    - `run_eqtm_with_lmem`:  eQTMs calculation using lmer
    - `run_eqtm_with_lmem_with_snps`: eQTMs calculation using lm and adjusting for meQTL SNP
- `03_eqtm_analysis`

### Results

The results are stored on the MPIP computational cluster: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/eqtms`
