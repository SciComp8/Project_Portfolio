suppressPackageStartupMessages(easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "foreach"))
source("./2022Dec8Version_ALiu/R/BMAseq.multi.postprob.norm.R")
source("./2022Dec8Version_ALiu/R/BMAseq.multi.postprob.MSout.norm.R")
source("./2022Dec8Version_ALiu/R/Bayesfactor.R")

var.pool0 <- c("BMI", "SEX") 
MS.default.train0 <- Modelspace(dat.pheno = dat.pheno.train, 
                                var.pool = var.pool0, 
                                max.nvar = 2, 
                                interaction = "BMI&SEX")
myModelSpace <- c(MS.default.train0$model.space,
                  paste0(myModelSpace0, "+AGE"),
                  paste0(myModelSpace0, "+MHABNWBC"),
                  paste0(myModelSpace0, "+AGE+MHABNWBC"))

MS.default.train <- Modelspace(dat.pheno = dat.pheno.train, 
                               var.pool = var.pool, # var.pool should not contain "BMI&SEX"
                               max.nvar = 4, # Maximum number of main-effect variables 
                               interaction = "BMI&SEX")

multiInt_ana_TMM_top_v2 <- function(seed.num = 999, threshold = 2000) {
  
  ##------set variables of interest------
  var.pool <- c("BMI", "AGE", "SEX", "MHABNWBC")  
  var.pool.add <- c(var.pool, "BMIxSEX")
  
  
  ##------Expanded BMAseq------
  output.multi.int1.train <- BMAseq.multi.postprob.MSout.norm(
    dat.expr.counts = dat.expr.train, 
    dat.pheno = MS.default.train$dat.pheno.new,
    model.space = myModelSpace, 
    cut.BF = 1)
  
  output.multi.int2.train <- BMAseq.multi.DEG(
    dat.pheno = output.multi.int1.train$dat.pheno.new,
    model.space = output.multi.int1.train$model.space, 
    post.modelprob = output.multi.int1.train$post.modelprob, 
    var.pool = var.pool,
    interact = T, 
    cut.FDR = 0.05)
  
  output.multi.int1.test <- BMAseq.multi.postprob.MSout.norm(
    dat.expr.counts = dat.expr.test, 
    dat.pheno = MS.default.test$dat.pheno.new,
    model.space = myModelSpace, 
    cut.BF = 1)
  
  output.multi.int2.test <- BMAseq.multi.DEG(
    dat.pheno = output.multi.int1.test$dat.pheno.new,
    model.space = output.multi.int1.test$model.space, 
    post.modelprob = output.multi.int1.test$post.modelprob, 
    var.pool = var.pool,
    interact = T, 
    cut.FDR = 0.05)
  
  BMAseq.eFDR.Main.train <- mclapply(1:length(var.pool),
                                     function(i) output.multi.int2.train$eFDR.Main[, i][order(output.multi.int2.train$eFDR.Main[, i])[1:threshold]],
                                     mc.cores = 4)
  
  BMAseq.eFDR.Main.test <- mclapply(1:length(var.pool),
                                         function(i) output.multi.int2.test$eFDR.Main[, i][order(output.multi.int2.test$eFDR.Main[, i])[1:threshold]],
                                         mc.cores = 4)
  
  names(BMAseq.eFDR.Main.train) = names(BMAseq.eFDR.Main.test) = var.pool
  
  BMAseq.eFDR.Interaction.train <- mclapply(1:length(var.pool),
                                            function(i) output.multi.int2.train$eFDR.Interaction[, i][order(output.multi.int2.train$eFDR.Interaction[, i])[1:threshold]],
                                            mc.cores = 4)
  
  BMAseq.eFDR.Interaction.test <- mclapply(1:length(var.pool),
                                           function(i) output.multi.int2.test$eFDR.Interaction[, i][order(output.multi.int2.test$eFDR.Interaction[, i])[1:threshold]],
                                           mc.cores = 4)
  
  names(BMAseq.eFDR.Interaction.train) = names(BMAseq.eFDR.Interaction.test) = var.pool
  
  save(dat.expr.train, dat.expr.test, dat.pheno.train, dat.pheno.test, var.pool, var.pool.add, 
       output.multi.int1.train, output.multi.int2.train, output.multi.int1.test, output.multi.int2.test, 
       BMAseq.eFDR.Main.train, BMAseq.eFDR.Main.test, BMAseq.eFDR.Interaction.train, BMAseq.eFDR.Interaction.test, 
       file = sprintf("../ApplicationData/derived/RandomSeed/Top%s/BMAseqMultiInt%s.RData", threshold, seed.num)) # Top%s needs to be created in advance; dir.create
}
