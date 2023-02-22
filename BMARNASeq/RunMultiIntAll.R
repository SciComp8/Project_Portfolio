suppressPackageStartupMessages(easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "foreach"))
source("./2022Dec8Version_ALiu/R/BMAseq.multi.postprob.norm.R")
source("./2022Dec8Version_ALiu/R/Bayesfactor.R")

multiInt_ana_TMM_top <- function(seed.num = 999, threshold = 2000) {
  
  ##------set variables of interest------
  var.pool <- c("BMI", "AGE", "SEX", "MHABNWBC")  
  var.pool.add <- c(var.pool, "BMIxSEX")
  
  
  ##------BMAseq------
  output.multi.int1.train <- BMAseq.multi.postprob.norm(
    dat.expr.counts = dat.expr.train,  
    dat.pheno = dat.pheno.train, 
    var.pool = var.pool, 
    max.nvar = 4, 
    interaction = "BMI&SEX", 
    cut.BF = 1)
  
  output.multi.int2.train <- BMAseq.multi.DEG(
    dat.pheno = output.multi.int1.train$dat.pheno.new,
    model.space = output.multi.int1.train$model.space, 
    post.modelprob = output.multi.int1.train$post.modelprob, 
    var.pool = var.pool,
    interact = T, 
    cut.FDR = 0.05)
  
  output.multi.int1.test <- BMAseq.multi.postprob.norm(
    dat.expr.counts = dat.expr.test, 
    dat.pheno = dat.pheno.test, 
    var.pool = var.pool, 
    max.nvar = 4, 
    interaction = "BMI&SEX", 
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
  
  names(BMAseq.eFDR.Main.train) = names(BMAseq.eFDR.Main.2000.test) = var.pool
  
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
       file = sprintf("../ApplicationData/derived/RandomSeed/Top%s/BMAseqMultiInt%s.RData", threshold, seed.num))
  
  
  ##------voom + limma------
  design.train <- model.matrix(~BMI + AGE + SEX + MHABNWBC + BMI*SEX, data = dat.pheno.train)
  design.test <- model.matrix(~BMI + AGE + SEX + MHABNWBC + BMI*SEX, data = dat.pheno.test)
  
  voom.train <- voom(counts = dat.expr.train,
                     design = design.train,
                     lib.size = colSums(dat.expr.train)*calcNormFactors(dat.expr.train)) 
  
  voom.fit.train <- lmFit(voom.train[["E"]],
                          design = design.train,
                          weights = voom.train[["weights"]])
  
  voom.eFDR.train <- mclapply(1:length(var.pool.add),
                              function(i) {
                                t <- voom.fit.train[["coefficients"]][, i + 1]/voom.fit.train[["stdev.unscaled"]][, i + 1]/voom.fit.train[["sigma"]]
                                p <- 2*pt(-abs(t), df = voom.fit.train[["df.residual"]])
                                return(qvalue(p)) },
                              mc.cores = 4L)
  
  voom.eFDR.train2 <- mclapply(1:length(var.pool.add),
                               function(i) {
                                 q.val <- voom.eFDR.train[[i]][["qvalues"]]
                                 return(q.val[order(q.val)[1:threshold]])
                               },
                               mc.cores = 4L)
  
  names(voom.eFDR.train) = names(voom.eFDR.train2) = var.pool.add 
  
  voom.test <- voom(counts = dat.expr.test,
                    design = design.test,
                    lib.size = colSums(dat.expr.test)*calcNormFactors(dat.expr.test))
  
  voom.fit.test <- lmFit(voom.test[["E"]],
                         design = design.test,
                         weights = voom.test[["weights"]])
  
  voom.eFDR.test <- mclapply(1:length(var.pool.add),
                             function(i) {
                               t <- voom.fit.test[["coefficients"]][, i + 1]/voom.fit.test[["stdev.unscaled"]][, i + 1]/voom.fit.test[["sigma"]]
                               p <- 2*pt(-abs(t), df = voom.fit.test[["df.residual"]])
                               return(qvalue(p)) },
                             mc.cores = 4L)
  
  voom.eFDR.test2 <- mclapply(1:length(var.pool.add),
                              function(i) {
                                q.val <- voom.eFDR.test[[i]][["qvalues"]]
                                return(q.val[order(q.val)[1:threshold]])
                              },
                                mc.cores = 4L)
  
  names(voom.eFDR.test) = names(voom.eFDR.test2) = var.pool.add 
  
  save(dat.expr.train, dat.expr.test, dat.pheno.train, dat.pheno.test, var.pool, var.pool.add, 
       design.train, design.test, voom.train, voom.fit.train, voom.eFDR.train, voom.eFDR.train2, 
       voom.test, voom.fit.test, voom.eFDR.test, voom.eFDR.test2, 
       file = sprintf("../ApplicationData/derived/RandomSeed/Top%s/VoomLimmaMultiInt%s.RData", threshold, seed.num))
  
  
  ##------voom + limma + eBayes------
  voom.train <- voom(counts = dat.expr.train,
                     design = design.train,
                     lib.size = colSums(dat.expr.train)*calcNormFactors(dat.expr.train))
  
  eBayes.fit.train <- lmFit(voom.train[["E"]],
                            design = design.train,
                            weights = voom.train[["weights"]]) %>% 
    eBayes()
  
  eBayes.eFDR.train <- mclapply(1:length(var.pool.add),
                                function(i) eBayes.fit.train[["p.value"]][, i + 1] %>% 
                                  qvalue(),
                                mc.cores = 4L)
  
  eBayes.eFDR.train2 <- mclapply(1:length(var.pool.add),
                                 function(i) {
                                   q.val <- eBayes.eFDR.train[[i]][["qvalues"]]
                                   return(q.val[order(q.val)[1:threshold]])
                                 },
                                 mc.cores = 4L)
  
  names(eBayes.eFDR.train) = names(eBayes.eFDR.train2) = var.pool.add 
  
  voom.test <- voom(counts = dat.expr.test,
                    design = design.test,
                    lib.size = colSums(dat.expr.test)*calcNormFactors(dat.expr.test))
  
  eBayes.fit.test <- lmFit(voom.test[["E"]],
                           design = design.test,
                           weights = voom.test[["weights"]]) %>% 
    eBayes()
  
  eBayes.eFDR.test <- mclapply(1:length(var.pool.add),
                               function(i) eBayes.fit.test[["p.value"]][, i + 1] %>% 
                                 qvalue(),
                               mc.cores = 4L)
  
  eBayes.eFDR.test2 <- mclapply(1:length(var.pool.add),
                                function(i) {
                                  q.val <- eBayes.eFDR.test[[i]][["qvalues"]]
                                  return(q.val[order(q.val)[1:threshold]]) 
                                },
                                mc.cores = 4L)
  
  names(eBayes.eFDR.test) = names(eBayes.eFDR.test2) = var.pool.add 
  
  save(dat.expr.train, dat.expr.test, dat.pheno.train, dat.pheno.test, var.pool, var.pool.add, 
       design.train, design.test, voom.train, eBayes.fit.train, eBayes.eFDR.train, eBayes.eFDR.train2, 
       voom.test, eBayes.fit.test, eBayes.eFDR.test, eBayes.eFDR.test2, 
       file = sprintf("../ApplicationData/derived/RandomSeed/Top%s/eBayesMultiInt%s.RData", threshold, seed.num))
  
  
  ##------edgeR------
  y.train <- DGEList(counts = dat.expr.train, 
                     lib.size = colSums(dat.expr.train)) %>% 
    calcNormFactors() %>% 
    estimateGLMTrendedDisp(design.train)
  
  edgeR.fit.train <- glmQLFit(y.train, design.train)
  
  qlf.train <- mclapply(1:length(var.pool.add),
                        function(i) glmQLFTest(edgeR.fit.train, coef = i + 1),
                        mc.cores = 4L)
  
  edgeR.eFDR.train <- mclapply(1:length(var.pool.add),
                               function(i) qlf.train[[i]][["table"]][["PValue"]] %>%
                                 qvalue(),
                               mc.cores = 4L)
  
  edgeR.eFDR.train2 <- mclapply(1:length(var.pool.add),
                                function(i) {
                                  q.val <- edgeR.eFDR.train[[i]][["qvalues"]]
                                  return(q.val[order(q.val)[1:threshold]]) 
                                },
                                mc.cores = 4L)
  
  edgeR.eFDR.GeneName.train <- mclapply(1:length(var.pool.add),
                                        function(i) {
                                          q.val <- edgeR.eFDR.train[[i]][["qvalues"]]
                                          return(rownames(dat.expr.train)[order(q.val)[1:threshold]])
                                        },
                                        mc.cores = 4L)
  
  names(qlf.train) = names(edgeR.eFDR.train) = names(edgeR.eFDR.train2) = names(edgeR.eFDR.GeneName.train) = var.pool.add 
  
  y.test <- DGEList(counts = dat.expr.test, 
                    lib.size = colSums(dat.expr.test)) %>% 
    calcNormFactors() %>% 
    estimateGLMTrendedDisp(design.test)
  
  edgeR.fit.test <- glmQLFit(y.test, design.test)
  
  qlf.test <- mclapply(1:length(var.pool.add),
                       function(i) glmQLFTest(edgeR.fit.test, coef = i + 1),
                       mc.cores = 4L)
  
  edgeR.eFDR.test <- mclapply(1:length(var.pool.add),
                              function(i) qlf.test[[i]][["table"]][["PValue"]] %>% 
                                qvalue(),
                              mc.cores = 4L)
  
  edgeR.eFDR.test2 <- mclapply(1:length(var.pool.add),
                               function(i) {
                                 q.val <- edgeR.eFDR.test[[i]][["qvalues"]]
                                 return(q.val[order(q.val)[1:threshold]]) 
                               },
                               mc.cores = 4L)
  
  edgeR.eFDR.GeneName.test <- mclapply(1:length(var.pool.add),
                                       function(i) {
                                         q.val <- edgeR.eFDR.test[[i]][["qvalues"]]
                                         return(rownames(dat.expr.test)[order(q.val)[1:threshold]])
                                       },
                                       mc.cores = 4L)
  
  names(qlf.test) = names(edgeR.eFDR.test) = names(edgeR.eFDR.test2) = names(edgeR.eFDR.GeneName.test) = var.pool.add 
  
  save(dat.expr.train, dat.expr.test, dat.pheno.train, dat.pheno.test, var.pool, var.pool.add, design.train, design.test, 
       y.train, edgeR.fit.train, edgeR.eFDR.train, edgeR.eFDR.train2, edgeR.eFDR.GeneName.train, 
       y.test, edgeR.fit.test, edgeR.eFDR.test, edgeR.eFDR.test2, edgeR.eFDR.GeneName.test, 
       file = sprintf("../ApplicationData/derived/RandomSeed/Top%s/edgeRMultiInt%s.RData", threshold, seed.num))
  
  
  ##------DESeq2------
  cts.train <- as.matrix(dat.expr.train)
  cts.test <- as.matrix(dat.expr.test)
  coldata.train <- dat.pheno.train
  coldata.test <- dat.pheno.test
  
  name.formula <- c("BMI_high_vs_low", "AGE_old_vs_young", "SEX_male_vs_female", 
                    "MHABNWBC_yes_vs_no", "BMIhigh.SEXmale")  
  
  lib.size <- colSums(cts.train)
  norm.factor <- calcNormFactors(cts.train, method = "TMM")
  size.factor <- lib.size*norm.factor/exp(mean(log(lib.size*norm.factor))) 
  
  dds <- DESeqDataSetFromMatrix(countData = cts.train, 
                                colData = coldata.train, 
                                design = ~BMI + AGE + SEX + MHABNWBC + BMI*SEX)
  sizeFactors(dds) <- size.factor
  
  res.train <- mclapply(1:length(var.pool.add), 
                        function(i) {
                          return(results(DESeq(dds), name = name.formula[i])) },
                        mc.cores = 4L)
  
  DESeq2.eFDR.train <- mclapply(1:length(var.pool.add),
                                function(i) qvalue(res.train[[i]][["pvalue"]])$qvalues,
                                mc.cores = 4L)
  
  DESeq2.eFDR.train2 <- mclapply(1:length(var.pool.add),
                                 function(i) {
                                   q.val <- DESeq2.eFDR.train[[i]]
                                   return(q.val[order(q.val)[1:threshold]])
                                 },
                                 mc.cores = 4L)
  
  DESeq2.eFDR.GeneName.train <- mclapply(1:length(var.pool.add),
                                         function(i) {
                                           q.val <- DESeq2.eFDR.train[[i]]
                                           return(rownames(cts.train)[order(q.val)[1:threshold]])
                                         },
                                         mc.cores = 4L)
  
  names(res.train) = names(DESeq2.eFDR.train) = names(DESeq2.eFDR.train2) = names(DESeq2.eFDR.GeneName.train) = var.pool.add 
  
  lib.size <- colSums(cts.test)
  norm.factor <- calcNormFactors(cts.test, method = "TMM")
  size.factor <- lib.size*norm.factor/exp(mean(log(lib.size*norm.factor))) 
  
  dds <- DESeqDataSetFromMatrix(countData = cts.test, 
                                colData = coldata.test, 
                                design = ~BMI + AGE + SEX + MHABNWBC + BMI*SEX)
  sizeFactors(dds) <- size.factor
  
  res.test <- mclapply(1:length(var.pool.add), 
                       function(i) {
                         return(results(DESeq(dds), name = name.formula[i])) }, 
                       mc.cores = 4L)
  
  DESeq2.eFDR.test <- mclapply(1:length(var.pool.add),
                               function(i) qvalue(res.test[[i]][["pvalue"]])$qvalues,
                               mc.cores = 4L)
  
  DESeq2.eFDR.test2 <- mclapply(1:length(var.pool.add),
                                function(i) {
                                  q.val <- DESeq2.eFDR.test[[i]]
                                  return(q.val[order(q.val)[1:threshold]])
                                },
                                mc.cores = 4L)
  
  DESeq2.eFDR.GeneName.test <- mclapply(1:length(var.pool.add),
                                        function(i) {
                                          q.val <- DESeq2.eFDR.test[[i]]
                                          return(rownames(cts.test)[order(q.val)[1:threshold]])
                                        },
                                        mc.cores = 4L)
  
  names(res.test) = names(DESeq2.eFDR.test) = names(DESeq2.eFDR.test2) = names(DESeq2.eFDR.GeneName.test) = var.pool.add
  
  save(cts.train, cts.test, coldata.train, coldata.test, var.pool, var.pool.add, 
       res.train, DESeq2.eFDR.train, DESeq2.eFDR.train2, DESeq2.eFDR.GeneName.train, 
       res.test, DESeq2.eFDR.test, DESeq2.eFDR.test2, DESeq2.eFDR.GeneName.test, 
       file = sprintf("../ApplicationData/derived/RandomSeed/Top%s/DESeq2MultiInt%s.RData", threshold, seed.num))
  
}
