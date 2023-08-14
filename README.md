# Molecular & clinical project-oriented programming

This repository compiles my programming work affiliated with collaborative and individual research projects, which span from statistical omics to healthcare research. 

## Table of contents

* [Bayesian model averaging for RNA-Seq count data](#Bayesian-model-averaging-for-RNA-Seq-count-data)
* [Statistical analysis for breast cancer surgery data](#Statistical-analysis-for-breast-cancer-surgery-data)
* [Differential expression analysis pipeline for single cell data](#Differential-expression-analysis-pipeline-for-single-cell-data)
* [Bayesian integrative clustering for multi-omics data](#Bayesian-integrative-clustering-for-multi-omics-data)
* [Lipidomics and transcriptomics analysis](#Lipidomics-and-transcriptomics-analysis)

## **Bayesian model averaging for RNA-Seq count data**
One of the fundamental research problems in bulk or single-cell transcriptomic studies is the identification of reliable differentially expressed genes (DEGs) between distinct sample groups (e.g., healthy and disease). Accurately identifying DEGs can offer reliable insights into biological processes underlying specific conditions, and further aid in developing diagnostic biomarkers and targeted therapeutics. This project evaluates the performance of our novel multi-model `BMAseq` and other existing single-model approaches: `DESeq2`, `edgeR`, and `limma` on inferring replicable DEGs under the univariable, multivaribles, and multivaribles with interaction terms model designs. The RNA-seq count and phenotypic data from GTEx are randomly divided into 50% training set and 50% test set. For evaluation, we focus on measures that are of practical interest to researchers when making decisions about selecting the appropriate model design to analyze their datasets. We investigate 1) number and rate of replicable DEGs at different ranking thresholds detected with different approaches; 2) consistency of replicable DEG detections between different approaches; 3) uniqueness of replicable DEG detections between different approaches; 4) major consistent detection ways of replicable DEGs; 5) biological insights of replicable DEGs that are uniquely and consistently identified by different approaches.

- [AnalysisPipeline.Rmd](BMARNASeq/AnalysisPipeline.Rmd)
- [GenerateCheck_TrainingTestData.Rmd](BMARNASeq/GenerateCheck_TrainingTestData.Rmd)

## **Statistical analysis for breast cancer surgery data**
**Multidisciplinary breast program in New York City** - This project uses the descriptive statistics, hypothesis tests, and multivariable regression modeling to delve deepper into the distributions of diverse factors such as demographics of the communities, insurance type, mammographic density, T and N stages, and tumor subtype in two hospitals, and to detect the predictors related to presenting at an early stage of breast cancer. By examining both aggregated and disaggregated racial/ethnic and hospital site data, we are interested in identifying the disparities in breast cancer diagnosis and treatment.
- [MDC_v1.Rmd](MDC/MDC_v1.Rmd)
- [MDC_v2.Rmd](MDC/MDC_v2.Rmd)
- [MDC_PredictiveModeling.Rmd](MDC/MDC_PredictiveModeling.Rmd)

**Triple negative breast cancer** - In this project, we evaluate the incidence and phenotypic distribution of subsequent ipsilateral or new primary contralateral breast cancer in patients with a history of triple negative breast cancer. The incidence rates of local recurrence and new primary breast cancer are estimated using the Kaplan-Meier approach. The cumulative incidence rates of local recurrence and new primary breast cancer with ER, PR, HER2 marker values are estimated using the cumulative incidence function, which accounts for the positive marker and the negative marker as competing risks. A by-product of this project is a pipeline that predicts the subsequent breast cancer event in 5 years using the machine learning approaches (e.g., random forest).

- [TNBC_v1.Rmd](TNBC/TNBC_v1.Rmd)
- [TNBC_v2.Rmd](TNBC/TNBC_v2.Rmd)
- [MachineLearningPipeline.py](TNBC/MachineLearningPipeline.py)
- [TNBC_StatisticalPlan.md](TNBC/TNBC_StatisticalPlan.md)

**Micrometastatic axillary nodal disease** - Most analytical work centers on the descriptive analysis, for example, estimate the number of patients undergoing completion axillary lymph node dissection (ALND), number of cases with additional positive nodes if completion ALND, and frequency of axillary recurrence stratified by whether ALND was performed and by whether radiation was delivered. The univariate survival analysis using the Cox-proportional hazards model is performed to interrogate the factors associated with the recurrence-free survival. Specially, recurrence-free survival time is defined as a composite endpoint which combines the time to recurrence, time to death, and time to last follow-up into a univariate time to the first event or censoring. 

- [Micromets.Rmd](Micromets.Rmd)

#### Note: The coding work related to other clinical projects can be found [here](https://github.com/anniliu7/Projectcoding/tree/main/Otherclinicalproject).
  
## **Differential expression analysis pipeline for single cell data**
The bioinformatics pipeline using the next-generation single cell data from the mouse B lymphocytes is developed to streamline the fundamental procedures of differentially expressed genes detection in single cell RNA-seq data. This includes the initial trimming of reads from raw FASTQ files using `Trim Galore`, quantifying the abundance of transcripts using `kallisto`, and performing the differential gene expression analysis using `DESeq2`. The pipeline is designed to be efficient in the command-line interface, enabling researchers to analyze and interpret large datasets of single cell gene expression data with the desirable speed, accuracy, and reproducibility.
  - [NGS_DEPipeline.Rmd](NGSBcell/NGS_DEPipeline.Rmd)
  - [NGS_TrimCount.sh](NGSBcell/NGS_TrimCount.sh)
  - [NGS_DEReport.Rmd](NGSBcell/NGS_DEReport.Rmd)

## **Bayesian integrative clustering for multi-omics data**
This project characterizes the shared patterns underlying the transcriptomic and lipidomic profiles in women at high risk or diagnosed with breast cancers and explores the relationships between these shared patterns and clinical and biological phenotypes. A novel Bayesian latent variable approach which jointly models the non-tumorous breast tissue transcriptomic expression levels and plasma lipidomic expression levels is constructed and further tuned with key parameters to identify clinically and biologically relevant clusters co-driven by transcriptomic and lipidomic features. These clusters and features are hidden by the traditional principal component analysis workflow.
  - [BayesClustering.Rmd](BayesClustering.Rmd)
 
## **Lipidomics and transcriptomics analysis**
The goal of this lipidomics and transcriptomics analysis is to test the hypothesis that lipid and gene molecules may be co-up- or co-down-regulated in different phenotypes (e.g., body mass index, trunk fat percent, adipocyte diameter states) and these differences may enable useful insights into the pathogenic mechanisms of breast cancer, to pave the way for developing more targeted prevention and intervention strategies for breast cancer patients. The main theme of analysis starts from processing lipidomic and transcriptomic expression data, analyzing the co-expression patterns in different phenotypes, to creating heatmap presentations. 
  - [LipidomicTranscriptomic.Rmd](LipidomicTranscriptomic.Rmd)

