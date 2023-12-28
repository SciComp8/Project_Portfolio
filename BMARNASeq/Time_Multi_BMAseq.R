#----Attach libraries and functions----  
library(microbenchmark)
library(benchmarkme)
library(parallel)
library(qvalue)
library(edgeR)
library(DESeq2)
library(BMAseq) 
source("BMAseq.multi.postprob.norm.R")
Bayesfactor <- BMAseq:::Bayesfactor
source("benchmark_model_test.R")

#----Load the benchmarking data----  
seed.i <- 8809678
dat.expr.train <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i)) 
dat.expr.test <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i)) 
dat.pheno.train <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i)) 
dat.pheno.test <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i)) 
vars.pool <- c("BMI", "AGE", "SEX", "MHABNWBC")  
name.formula <- c("BMI_high_vs_low", "AGE_old_vs_young", "SEX_male_vs_female", "MHABNWBC_yes_vs_no") 
design.train <- model.matrix(~BMI + AGE + SEX + MHABNWBC, data = dat.pheno.train)
design.test <- model.matrix(~BMI + AGE + SEX + MHABNWBC, data = dat.pheno.test)

#----View machine's specifications----  
get_cpu()$model_name 
get_cpu()$no_of_cores 
get_ram() 
get_platform_info()$OS.type
get_r_version()$platform 
get_r_version()$version.string 

#----Time multi-model BMAseq, multivariable DESeq2, multivariable edgeR, multivariable eBayes, multivariable voom.limma----
microbenchmark(multi_TMM_top_BMAseq(seed.val = seed.i, use.train = TRUE), times = 10)

microbenchmark(multi_TMM_top_DESeq2(seed.val = seed.i, use.train = TRUE), times = 10)

microbenchmark(multi_TMM_top_edgeR(seed.val = seed.i, use.train = TRUE), times = 10)

microbenchmark(multi_TMM_top_eBayes(seed.val = seed.i, use.train = TRUE), times = 10)

microbenchmark(multi_TMM_top_voom.limma(seed.val = seed.i, use.train = TRUE), times = 10)
