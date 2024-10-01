## **psychoGxE:** joint impact of environmental and genetic determinants on the aetiology of stress-related psychiatric disorders

- [Objectives](#objectives)
- [Data](#data)
- [Analysis](#analysis)
  - [1. DNA methhylation (DNAm) profiling](#1-dna-methhylation-dnam-profiling)
  - [2. Smoking score estimation](#2-smoking-score-estimation)
  - [3. Differential DNA methylation analysis](#3-differential-dna-methylation-analysis)
  - [4. Methylation quantitative trait loci (meQTL) analysis](#4-methylation-quantitative-trait-loci-meqtl-analysis)
  - [5. Expression quantitative trait methylation (eQTM) analysis](#5-expression-quantitative-trait-methylation-eqtm-analysis)
  - [6. Functional genomic annotation of differential DNAm](#6-functional-genomic-annotation-of-differential-dnam)
  - [7. Cell type composition estimation](#7-cell-type-composition-estimation)
  - [8. Cell-type specific analysis](#8-cell-type-specific-analysis)
  - [9. Functional enrichment analysis](#9-functional-enrichment-analysis)
  - [10. Histone mark enrichment analysis](#10-histone-mark-enrichment-analysis)
  - [11. Genome-wide association studies (GWAS) enrichment analysis](#11-genome-wide-association-studies-gwas-enrichment-analysis)
  - [12. Phenome-wide association studies (PheWAS) analysis](#12-phenome-wide-association-studies-phewas-analysis)
- [Results and additional information](#results-and-additional-information)
- [Authors](#Authors)

## Objectives

The main objectives of the study are to explore the genetic variants that moderate stress-induced environmental patterns in an immune cell-type-specific manner; and  assess how these epigenetic, i.e. environmental, and genetic biomarkers relate to the risk for psychiatric disorders.

## Data

Participants comprised **289** Caucasian individuals of the Max Planck Institute of Psychiatry (MPIP) in Munich. Of the participants, 

+ 93 women and 196 men
+ 129 (81 men, 48 women) were treated for major dipressive disorders (MDD)
+ 160 (115 men, 45 women) were healthy controls with no history of a depressive disorder. 

Baseline whole blood samples were obtained at 6 pm after two hours of fasting and abstention from coffee and physical activity. Subjects then received 1.5 mg oral dexamethasone, and a second blood draw was performed at 9 pm, three hours after dexamethasone ingestion.

The available multiomics data for the currect project:

- Methylation (DNAm)
- Gene-expression
- Genotype
- Phenotype 

Data are stored on the MPIP computational cluster:

## Analysis

The following figure displays a flow of the data collection and statistical analyses outlined below:

![workflow](https://github.com/ahryho/psychoGxE/blob/master/materials/figures/dex-stim-workflow.jpg)

The analysis includes the follwong steps:

### 1. DNA methhylation (DNAm) profiling
   
For the methodology and scripts, please refer to [DNAm QC](https://github.com/ahryho/dex-stim-human-dna-methyl-qc#dex-stimulated-dnam-arrays-processing).

### 2. Smoking score estimation
   
For the methodology and scripts, please refer to [smoking score calculation](https://github.com/ahryho/dex-stim-human-smoking-score).

### 3. Differential DNA methylation analysis

For the methodology and scripts, please refer to [dDNAm](https://github.com/ahryho/psychogxe/tree/master/code/methylation/).

### 4. Methylation quantitative trait loci (meQTL) analysis

For the methodology and scripts, please refer to the [meqtl folder](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/meqtl).

### 5. Expression quantitative trait methylation (eQTM) analysis

For the methodology and scripts, please refer to [eQTM analysis](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/eqtm).

### 6. Functional genomic annotation of differential DNAm

Differential DNA methylated positions were mapped to their genomic and gene location according to Illumina’s annotation using the R package minfi. To assess the distribution of CpG sites across the genome, the number of probes was calculated and weighted by the corresponding chromosome length for each chromosome. To assess the relation of CpGs to the gene and genomic regions, the distance from all CpGs to each of their annotated genes was calculated, and the closest genes were taken. In the case of multiple regions mapped to a single CpG, the CpG was added to each region.

For more details, please refer to [scripts](https://github.com/ahryho/psychoGxE/tree/master/code/methylation/02_dma/01_dmp).

### 7. Cell type composition estimation

The blood cell-type components were predicted for DNA methylation data on whole blood samples from the MPIP cohort. The model proposed by Salas et al. and R package EpiDISH were used to estimate the proportions of six main cell types in whole blood (CD4+ T cells, CD8+ T cells, monocytes, B cells, granulocytes, and natural killer cells) as well as subtypes of T and B cells (naïve, memory, and regulatory CD4+ T cells as well as naïve and memory CD8+ T cells and naïve and memory B cells).

For more details, please refer to [methodology](https://github.com/ahryho/dex-stim-human-dna-methyl-qc#9-cell-types-estimation) and [scripts](https://github.com/ahryho/dex-stim-human-dna-methyl-qc/tree/master/09_estimate_cell_proportion).

### 8. Cell-type specific analysis

For the methodology and scripts, please refer to the [cell-type specific analysis](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/cell_type_specific_analysis).

### 9. Functional enrichment analysis

For the methodology and scripts, please refer to the [functional enrichment of meQTL CpGs](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_cpg_annotation#functional-enrichment-analysis).

### 10. Histone mark enrichment analysis

For the methodology and scripts, please refer to the histone mark enrichment analysis of [meQTL CpGs](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_cpg_annotation#histone-mark-enrichment-analysis) and [meQTL SNPs](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation#histone-mark-enrichment-analysis).

### 11. Genome-wide association studies (GWAS) enrichment analysis

For the methodology and scripts, please refer to the [GWAS enrichment analysis](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/meqtl/04_me-qtl_analysis/04_04_me-qtl_snp_annotation#gwas-enrichment-analysis).

### 12. Phenome-wide association studies (PheWAS) analysis

For the methodology and scripts, please refer to the [PheWAS analysis](https://github.com/ahryho/psychoGxE/tree/master/code/integrative/phewas).

## Results and additional information

Project outcomes are stored on the MPIP computational cluster.

For additional information on methodology, results and limitations, please contact the owner of this repository.

## Authors
[Anastasiia Hryhorzhevska](https://www.linkedin.com/in/ahryhorzhevska)
