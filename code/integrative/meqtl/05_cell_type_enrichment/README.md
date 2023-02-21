# Cell-type specific analysis
The cell-type specific analysis was conducted using the linear model on the baseline and GR-response (post-dexamethasone) DNA methylation beta values. To assess possible violation of regression assumptions56, variance inflation factors (VIFs) were calculated using the R package car57 for all coefficients in regression equations for each CpG. VIFs measure how much the variances of the regression coefficients are inflated compared to no correlation between predictor variables. The percentage of variance explained by the other independent variables can be calculated using the formula (1-1/VIF)×100. A VIF of 5 or less is generally considered acceptable56. Regression models were re-run, excluding the cell type with the highest VIF (Neutrophils) to reduce multicollinearity. In addition, each model was adjusted for covariates: sex, age, BMI, depression status (MDD), smoking score and first two PCs from genotypic data:

CpG ∼ β_0 + β_{1-11} × BCC_{1-11} + β_12 × Sex + β_13 × Age + β_14 × BMI + β_15 × MDD + β_{16,17} × PC_{1,2} + β_18 × SS + ϵ 

### Assessment of the relationship between CpG sites and Neutrophils

To assess the relationship between CpG sites and Neutrophils, linear regression models were built, excluding those cell types highly correlated with Neutrophils. The correlation between blood cell counts was measured using the Pearson correlation coefficient. 

The FDR threshold of 0.05 was used to identify the number of significant baseline and GR-response CpG-blood cell type associations. 

### Assessment whether the meQTL CpGs are cell-type-specific

To assess whether the meQTL CpGs are cell-type-specific, the significant baseline and GR-response CpG sites among those in common were compared. In addition, the Kolmogorov-Smirnov test was used to evaluate the significance between GR-response and baseline for each blood cell type.
