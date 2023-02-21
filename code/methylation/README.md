## Differential DNA methylation analysis

To observe dexamethasone-dependent changes in DNA methylation, a linear mixed effect model was applied to the filtered, normalised, and batch-corrected methylation beta values controlling for sex, age, BMI, depression status (MDD), first three SVs derived from DNA methylation data, and first two PCs which describes the genetic variability of the study cohort. For each CpG site, the following model was tested using the R package lme4: 

CpG ∼ β_0 + β_1 × Sex + β_2 × Age + β_3 × BMI + β_4 × MDD + β_(5-7)× SV_(5-7) + β_8,9 × PC1,2 + γ × DEX + ϵ

To identify the number of significant GR-response differentially methylated positions (DMPs), the false discovery rate, FDR < 0.01, and fold change, FC < 0.02, were used. Of note, the FC was calculated as the difference in DNA methylation between post-dexamethasone and baseline.
