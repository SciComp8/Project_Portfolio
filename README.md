# Coding collection of diverse research projects

Here is a summary of my coding work affiliated with multiple collaborative and individual research projects.

## **Bayesian model averaging for bulk RNA-Seq data (use R)**
  - The goal of this project is to evaluate the performance  of novel BMAseq and other three mainstream approaches: voom + limma, edgeR, and DESeq2 on identifying differentially expressed genes in the setting of univariable and multivaribles. All subjects are randomly divided into 50% training set and 50% test set. Each method of BMAseq,  voom + limma, edgeR, and DESeq2 is applied to construct the univariate and multivariate models for the following variables: BMI, AGE, SEX, MHABNWBC, MHARTHTS, and MHCVD. The approach-specific Venn diagram panel is plotted to show the proportion of common differentially expressed genes (cDEGs) associated with each variable between the training and test set. The proportion of cDEGs associated with each variable is compared across four approaches.

    - [UniBMA_voom_limma_edgeR_DESeq2.Rmd](UniBMA_voom_limma_edgeR_DESeq2.Rmd)
    - [MultiBMA_voom_limma_edgeR_DESeq2.Rmd](MultiBMA_voom_limma_edgeR_DESeq2.Rmd)
  
- **Differential expression analysis pipeline for single cell data (use Python-based snakefile and R)**
  - [NGS_DEPipeline.Rmd](NGS_DEPipeline.Rmd)

- **Bayesian integrative clustering for multi-omics data (use R)**
  - [BayesClustering.Rmd](BayesClustering.Rmd)
 
- **Lipidomics and transcriptomics analysis (use R)**
  - [LipidomicTranscriptomic.Rmd](LipidomicTranscriptomic.Rmd)

- **Statistical analysis for breast cancer surgery clinical data (use R and Python)**
  - [MDCAsian.Rmd](MDCAsian.Rmd)
  - [TNBC.Rmd](TNBC.Rmd): The goal of this project is to 
  - [Micromets.Rmd](Micromets.Rmd)
