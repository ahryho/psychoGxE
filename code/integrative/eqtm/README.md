## Expression quantitative trait methylation analysis

The baseline/GR-response expression quantitative trait methylation (eQTM) analysis was conducted on 11,944 baseline/GR-response gene expressions and baseline/GR-response 740,357 CpG sites using linear regression. Each model was adjusted for covariates: sex, age, BMI, depression status (MDD), first three SVs derived from DNA methylation data, first two PCs from genotypic data and Smoking score (SS):  

GEX $∼$ $β_0 + β_CpG \times$ CpG $+$ $β_1 \times$ Sex $+$ $β_2  Age + β_3 × BMI + β_4 × MDD + β_{5-7} × SV_{5-7} + β_8,9 × PC_{1,2} + β_10 × SS + ϵ 

The eQTM analysis aims to investigate whether the DNA methylation influences gene expression and the integration of the results with meQTL analysis. Therefore, the nominal p-value of 0.05 was used to identify the number of significant baseline and GR-response eQTM pairs.
