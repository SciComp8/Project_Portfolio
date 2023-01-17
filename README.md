# Research project-oriented programming

Here lists my coding work affiliated with multifaceted collaborative and individual research projects spanning from statistical omics to clinical research.


## **Bayesian model averaging for bulk RNA-Seq data (use R programming)**
The goal of this project is to evaluate the performance  of our novel BMAseq and other three existing biostatistical and bioinformatics approaches: voom + limma, edgeR, and DESeq2 on identifying differentially expressed genes in the setting of univariable, multivaribles, and multivaribles with interaction. All subjects are randomly divided into 50% training set and 50% test set. Each method of BMAseq,  voom + limma, edgeR, and DESeq2 is applied to construct the univariate, multivariate, multivariate with interaction models for the six phenotypes: BMI, AGE, SEX, MHABNWBC, MHARTHTS, and MHCVD. The approach-specific Venn diagram panel is plotted to intuitively reveal the absolute number and relative proportion of common differentially expressed genes (cDEGs) associated with each variable between the training and test set. The distribution of cDEGs associated with each variable is compared across four approaches.

- [UniBMA_voom_limma_edgeR_DESeq2.Rmd](UniBMA_voom_limma_edgeR_DESeq2.Rmd)
- [UniBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd](UniBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd)
- [MultiBMA_voom_limma_edgeR_DESeq2.Rmd](MultiBMA_voom_limma_edgeR_DESeq2.Rmd)
- [MultiBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd](MultiBMA_voom_limma_edgeR_DESeq2_Top2000.Rmd)
- [Interaction_multiBMA_voom_limma_edgeR_DESeq2.Rmd](Interaction_multiBMA_voom_limma_edgeR_DESeq2.Rmd)
  
## **Differential expression analysis pipeline for single cell data (use Python-based snakemake workflow and R programming)**
A bioinformatics pipeline using the next-generation single cell data from the mouse B lymphocytes is developed for the purpose of streamlining several steps involving trimming reads from raw FASTQ files, quantifying abundances of transcripts from RNA-Seq data, and performing downstream DESeq2-based differential gene expression analysis. 
  - [NGS_DEPipeline.Rmd](NGS_DEPipeline.Rmd)

## **Bayesian integrative clustering for multi-omics data (use R programming)**
This project characterizes the shared patterns underlying the transcriptomic and lipidomic profiles in women at high risk or diagnosed with breast cancers and explores the relationships between these shared patterns and clinical and biological phenotypes. A novel Bayesian latent variable approach which jointly models the non-tumorous breast tissue transcriptomic expression levels and plasma lipidomic expression levels is constructed and further tuned with key parameters to identify clinically and biologically relevant clusters co-driven by transcriptomic and lipidomic features.
  - [BayesClustering.Rmd](BayesClustering.Rmd)
 
## **Lipidomics and transcriptomics analysis (use R programming)**
The goal of this lipidomics and transcriptomics analysis is to test the hypothesis where lipid and gene molecules may be co-up- or co-down-regulated in different phenotypes (e.g., body mass index, trunk fat percent, adipocyte diameter states) and these differences may enable useful insights into the pathogenic mechanisms of breast cancer, ultimately for the purpose of developing more targeted prevention and intervention strategies for breast cancer patients. The analysis starts from raw lipidomic and transcriptomic data preprocessing, statistical analysis, and heatmap visualization. 
  - [LipidomicTranscriptomic.Rmd](LipidomicTranscriptomic.Rmd)

## **Statistical analysis for breast cancer surgery clinical data (use R/Python programming)**
The TNBC project aims to evaluate the incidence and phenotypic distribution of subsequent ipsilateral or new primary contralateral breast cancer in patients with a history of triple negative breast cancer. The incidence rates of local recurrence and new primary breast cancer are estimated in 1, 2, 3, 5, 10 years, separately, using the Kaplan-Meier approach. The cumulative incidence rates of local recurrence and new primary breast cancer with ER, PR, HER2 marker values are estimated in 1, 2, 3, 5, 10 years, separately, using the cumulative incidence function which considers the positive marker and the negative marker (e.g., ER positive vs ER negative) as competing risks. A by-product of this project is to construct a pipeline which predicts the subsequent breast cancer event in 5 years using the existing machine learning approaches (e.g. random forest).

- [TNBC.Rmd](TNBC.Rmd)
- [MachineLearningPipeline.py](MachineLearningPipeline.py)

In the micrometastatic axillary nodal disease project, most analytical work centers on the descriptive analysis, for example, estimate the number of patients undergoing completion axillary lymph node dissection (ALND), number of cases with additional positive nodes if completion ALND, and frequency of axillary recurrence stratified by whether ALND was performed and by whether radiation was delivered. The univariate survival analysis using the Cox-proportional hazards model is performed to interrogate the factors associated with the recurrence-free survival. Specially, recurrence-free survival time is defined as a composite endpoint which combines the time to recurrence, time to death, and time to last follow-up into a univariate time to the first event or censoring. 

- [Micromets.Rmd](Micromets.Rmd)

In the MDC Asian project, the breast cancer patterns such as frequency of screen-detected disease, clinical T and N stages, and tumor subtype are evaluated using the aggregated and disaggregated racial/ethnic and community demographic data at two hospitals.
- [MDCAsian.Rmd](MDCAsian.Rmd)
