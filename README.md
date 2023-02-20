## **psychoGE:** joint impact of environmental and genetic determinants on the aetiology of stress-related psychiatric disorders

The main objectives of the study are to explore the genetic variants that moderate stress-induced environmental patterns in an immune cell-type-specific manner; and  assess how these epigenetic, i.e. environmental, and genetic biomarkers relate to the risk for psychiatric disorders.

## Data

Participants comprised **289** Caucasian individuals of the Max Planck Institute of Psychiatry (MPIP) in Munich. Of the participants, 

+ 93 women and 196 men
+ 129 (81 men, 48 women) were treated for major dipressive disorders (MDD)
+ 160 (115 men, 45 women) were healthy controls with no history of a depressive disorder. 

Baseline whole blood samples were obtained at 6 pm after two hours of fasting and abstention from coffee and physical activity. Subjects then received 1.5 mg oral dexamethasone, and a second blood draw was performed at 9 pm, three hours after dexamethasone ingestion.

The available multiomics data for the currect project:

- Methylation
- Gene-expression
- Genotype
- Phenotype 

For detailed data overview, please refer to the ```/code/integrative/data_overview```

## Analysis

The analysis includes the follwong steps:

1. DNA methhylation (DNAm) profiling
   
   The workflow and scripts are located: [dnam qc](https://github.com/ahryho/dex-stim-human-dna-methyl-qc)

2. Smoking score estimation
3. Differential DNA methylation analysis
4. Methylation quantitative trait loci (meQTL) analysis
5. Expression quantitative trait methylation (eQTM) analysis
6. Functional genomic annotation of differential DNAm
7. Cell type composition estimation
8. Cell-type specific analysis
9.  Functional enrichment analysis
10. Histone mark enrichment analysis
11. Genome-wide association studies (GWAS) enrichment analysis
12. Phenome-wide association studies (PheWAS) analysis