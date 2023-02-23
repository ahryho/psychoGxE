## Differential DNA methylation analysis

To observe dexamethasone-dependent changes in DNA methylation, a linear mixed effect model (LMEM) was applied to the filtered, normalised, and batch-corrected methylation beta values controlling for sex, age, BMI, depression status (MDD), first three surrogate variabels (SVs) derived from DNA methylation data, and first two principal components (PCs) which describes the genetic variability of the study cohort. For each CpG site, the following model was tested using the R package lme4: 

CpG $∼ β_0+β_1 \times$ Sex $+$ $β_2 \times$ Age $+$ $β_3 \times$ BMI $+$ $β_4 \times$ MDD $+$ $β_{5-7} \times$ $SV_{5-7}$ $+$ $β_{8,9} \times$ $PC_{1,2}$ $+$ $γ \times$ DEX $+$ $ϵ$

To identify the number of significant GR-response differentially methylated positions (DMPs), the false discovery rate, **FDR $< 0.01$**, and fold change, **FC** $< 0.02$, were used. Of note, the FC was calculated as the difference in DNA methylation between post-dexamethasone and baseline.

### Folders structure

- `00_prepare the data`: the script for phenotype data to identify differentially methylated regions (DMRs)
- `01_lmem_dnam`: the scripts for DMPs identifications. 6 models were tested:
  
    - `lmem_bcc`: controlling for sex, age, BMI, MDD, blood cell counts (BCC) derived from DNAm
    - `lmem_bcc_cellcode_pcs`: controlling for sex, age, BMI, MDD, BCC derived from gene-expression, and genotype PCs
    - `lmem_bcc_pcs`: controlling for sex, age, BMI, MDD, BCCs derived from DNAm, genotype PCs
    - `lmem_no_bcc`: controlling for sex, age, BMI, MDD
    - `lmem_svs_pcs`: controlling for sex, age, BMI, MDD, SVs derived from DNAm, genotype PCs
    - `lmem_svs_pcs_smoke`: controlling for sex, age, BMI, MDD, SVs derived from DNAm, genotype PCs, and smoking score 
- `02-dma`: script for identification of DMPs and DMRs, and gene set enrichment analysis using rGREAT 
- `03_additional_statistics`: additional calculations

### Results

The results are stored on the MPIP computational cluster: `/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/methylation/`
