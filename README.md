# Coding collection of diverse research projects

Here is a summary of my coding work affiliated with multiple collaborative and individual research projects.

## **Bayesian model averaging for bulk RNA-Seq data (use R)**
  - The goal of this project is to evaluate the performance  of our novel BMAseq and other three existing biostatistical and bioinformatics approaches: voom + limma, edgeR, and DESeq2 on identifying differentially expressed genes in the setting of univariable, multivaribles, and multivaribles with interaction. All subjects are randomly divided into 50% training set and 50% test set. Each method of BMAseq,  voom + limma, edgeR, and DESeq2 is applied to construct the univariate, multivariate, multivariate with interaction models for the six phenotypes: BMI, AGE, SEX, MHABNWBC, MHARTHTS, and MHCVD. The approach-specific Venn diagram panel is plotted to intuitively reveal the proportion of common differentially expressed genes (cDEGs) associated with each variable between the training and test set. The proportion of cDEGs associated with each variable is compared across four approaches.

    - [UniBMA_voom_limma_edgeR_DESeq2.Rmd](UniBMA_voom_limma_edgeR_DESeq2.Rmd)
    - [MultiBMA_voom_limma_edgeR_DESeq2.Rmd](MultiBMA_voom_limma_edgeR_DESeq2.Rmd)
  
## **Differential expression analysis pipeline for single cell data (use Python-based snakefile and R)**
  - [NGS_DEPipeline.Rmd](NGS_DEPipeline.Rmd)

## **Bayesian integrative clustering for multi-omics data (use R)**
  - [BayesClustering.Rmd](BayesClustering.Rmd)
 
## **Lipidomics and transcriptomics analysis (use R)**
  - [LipidomicTranscriptomic.Rmd](LipidomicTranscriptomic.Rmd)

## **Statistical analysis for breast cancer surgery clinical data (use R and Python)**
  - The TNBC project aims to evaluate the cumulative incidence and phenotypic distribution of subsequent ipsilateral or new primary contralateral breast cancer in patients with a history of triple negative breast cancer. 
    - [TNBC.Rmd](TNBC.Rmd)
  - [MDCAsian.Rmd](MDCAsian.Rmd)
  - [Micromets.Rmd](Micromets.Rmd)
