suppressPackageStartupMessages(easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "foreach"))
source("./multiInt_ana_TMM_top_v2.R") 

threshold.vec <- seq(1000, 5000, 1000)
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
for (threshold.i in threshold.vec) {
  for (seed.i in seed.vec) {
    dat.expr.train <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i))
    dat.expr.test <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i))
    dat.pheno.train <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i))
    dat.pheno.test <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i))
    multiInt_ana_TMM_top_v2(seed.num = seed.i, threshold = threshold.i)
  }
}
