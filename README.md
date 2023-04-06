# Molecular & clinical project-oriented programming

This is a compiled list of my programming work affiliated with collaborative and individual research projects, which cover the topics from statistical omics to healthcare research. 

## Table of contents

* [Bayesian model averaging for RNA-Seq count data](#Bayesian-model-averaging-for-RNA-Seq-count-data)
* [Differential expression analysis pipeline for single cell data](#Differential-expression-analysis-pipeline-for-single-cell-data)
* [Bayesian integrative clustering for multi-omics data](#Bayesian-integrative-clustering-for-multi-omics-data)
* [Lipidomics and transcriptomics analysis](#Lipidomics-and-transcriptomics-analysis)
* [Statistical analysis for breast cancer surgery clinical data](#Statistical-analysis-for-breast-cancer-surgery-clinical-data)



## **Bayesian model averaging for RNA-Seq count data**
Accurately identifying differentially expressed genes(DEGs) can offer reliable insights into biological processes underlying specific conditions, and further aid in developing diagnostic biomarkers and targeted therapeutics. The goal of this project is to evaluate the performance of our novel `BMAseq` and other existing single-model approaches: `DESeq2`, `edgeR`, `voom + limma + eBayes`, and `voom + limma` on inferring replicable DEGs in the setting of univariable, multivaribles, and multivaribles with interaction terms. The RNA-seq count and phenotypic data from GTEx are randomly divided into 50% training set and 50% test set. Each approach is used to build the univariate, multivariate, multivariate with interaction models for each phenotype. The approach-specific Venn diagram panel and boxplot panel are plotted to reveal the number and rate of replicable DEGs associated with each phenotype between the training and test set. 

- [AnalysisPipeline.Rmd](BMARNASeq/AnalysisPipeline.Rmd)
- [UniBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd](BMARNASeq/UniBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd)
- [MultiBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd](BMARNASeq/MultiBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd)
- [Interaction_multiBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd](BMARNASeq/Interaction_multiBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd)
- [GenerateCheck_TrainingTestData.Rmd](BMARNASeq/GenerateCheck_TrainingTestData.Rmd)
  
## **Differential expression analysis pipeline for single cell data**
The bioinformatics pipeline developed for the next-generation single cell data from the mouse B lymphocytes aims to streamline the entire process of differential gene expression analysis for single cell RNA-seq data. This involves several key phases, including the initial trimming of reads from raw FASTQ files using `Trim Galore`, quantifying the abundance of transcripts using `kallisto`, and finally, performing the differential gene expression analysis using `DESeq2`. The pipeline is designed to be efficient in the command-line interface, enabling researchers to analyze and interpret large datasets of single cell gene expression data with the desirable speed, accuracy, and reproducibility.
  - [NGS_DEPipeline.Rmd](NGS_DEPipeline.Rmd)

## **Bayesian integrative clustering for multi-omics data**
This project characterizes the shared patterns underlying the transcriptomic and lipidomic profiles in women at high risk or diagnosed with breast cancers and explores the relationships between these shared patterns and clinical and biological phenotypes. A novel Bayesian latent variable approach which jointly models the non-tumorous breast tissue transcriptomic expression levels and plasma lipidomic expression levels is constructed and further tuned with key parameters to identify clinically and biologically relevant clusters co-driven by transcriptomic and lipidomic features.
  - [BayesClustering.Rmd](BayesClustering.Rmd)
 
## **Lipidomics and transcriptomics analysis**
The goal of this lipidomics and transcriptomics analysis is to test the hypothesis that lipid and gene molecules may be co-up- or co-down-regulated in different phenotypes (e.g., body mass index, trunk fat percent, adipocyte diameter states) and these differences may enable useful insights into the pathogenic mechanisms of breast cancer, to pave the way for developing more targeted prevention and intervention strategies for breast cancer patients. The main theme of analysis starts from processing lipidomic and transcriptomic expression data, analyzing the co-expression patterns in different phenotypes, to creating heatmap presentations. 
  - [LipidomicTranscriptomic.Rmd](LipidomicTranscriptomic.Rmd)

## **Statistical analysis for breast cancer surgery clinical data**
**TNBC project** - This project aims to evaluate the incidence and phenotypic distribution of subsequent ipsilateral or new primary contralateral breast cancer in patients with a history of triple negative breast cancer. The incidence rates of local recurrence and new primary breast cancer are estimated in 1, 2, 3, 5, 10 years, separately, using the Kaplan-Meier approach. The cumulative incidence rates of local recurrence and new primary breast cancer with ER, PR, HER2 marker values are estimated in 1, 2, 3, 5, 10 years, separately, using the cumulative incidence function which considers the positive marker and the negative marker (e.g., ER positive vs ER negative) as competing risks. A by-product of this project is to construct a pipeline which predicts the subsequent breast cancer event in 5 years using the existing machine learning approaches (e.g., random forest).

- [TNBC_v1.Rmd](TNBC/TNBC_v1.Rmd)
- [TNBC_v2.Rmd](TNBC/TNBC_v2.Rmd)
- [MachineLearningPipeline.py](TNBC/MachineLearningPipeline.py)
- [TNBC_StatisticalPlan.md](TNBC/TNBC_StatisticalPlan.md)

**Micrometastatic axillary nodal disease project** - In this project, most analytical work centers on the descriptive analysis, for example, estimate the number of patients undergoing completion axillary lymph node dissection (ALND), number of cases with additional positive nodes if completion ALND, and frequency of axillary recurrence stratified by whether ALND was performed and by whether radiation was delivered. The univariate survival analysis using the Cox-proportional hazards model is performed to interrogate the factors associated with the recurrence-free survival. Specially, recurrence-free survival time is defined as a composite endpoint which combines the time to recurrence, time to death, and time to last follow-up into a univariate time to the first event or censoring. 

- [Micromets.Rmd](Micromets.Rmd)

**MDC Asian project** - This project aims to delve deeper into breast cancer patterns by analyzing various factors such as disease frequency, clinical T and N stages, tumor subtype, and the demographics of the communities served by two hospitals. By examining both aggregated and disaggregated racial/ethnic and community demographic data, we identify the disparities in breast cancer diagnosis and treatment. 
- [MDCAsian.Rmd](MDCAsian.Rmd)
