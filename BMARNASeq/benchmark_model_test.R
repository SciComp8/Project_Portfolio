multi_TMM_top_BMAseq <- function(seed.val = 999, use.train = TRUE) {
  run_BMAseq <- function(dat.expr, dat.pheno, threshold.val) {
    output.multi.s1 <- BMAseq.multi.postprob.norm(dat.expr.counts = dat.expr,  
                                                  dat.pheno = dat.pheno, 
                                                  var.pool = vars.pool, 
                                                  max.nvar = 4, 
                                                  cut.BF = 1)
    output.multi.s2 <- BMAseq.multi.DEG(postprob.output = output.multi.s1, 
                                        cut.FDR = 1) 
    eFDR.Main <- mclapply(1:length(vars.pool),
                          function(i) output.multi.s2$eFDR[, i][order(output.multi.s2$eFDR[, i])[1:threshold.val]],
                          mc.cores = 10) 
    names(eFDR.Main) = vars.pool
    return(eFDR.Main)
  }
  
  BMAseq.eFDR.Main.train <- BMAseq.eFDR.Main.test <- NULL
  
  if (use.train) {
    BMAseq.eFDR.Main.train <- run_BMAseq(dat.expr.train, dat.pheno.train, 5000)
  } else {
    BMAseq.eFDR.Main.train <- run_BMAseq(dat.expr.train, dat.pheno.train, 5000)
    BMAseq.eFDR.Main.test <- run_BMAseq(dat.expr.test, dat.pheno.test, 5000)
  }
}

multi_TMM_top_DESeq2 <- function(seed.val = 999, use.train = TRUE) {
  run_DESeq2 <- function(cts, coldata, threshold.val) {
    lib.size <- colSums(cts)
    norm.factor <- calcNormFactors(cts, method = "TMM")
    size.factor <- lib.size * norm.factor / exp(mean(log(lib.size * norm.factor)))
    
    dds <- DESeqDataSetFromMatrix(countData = cts, 
                                  colData = coldata, 
                                  design = ~BMI + AGE + SEX + MHABNWBC)
    sizeFactors(dds) <- size.factor
    
    res <- mclapply(1:length(vars.pool), 
                    function(i) {
                      return(results(DESeq(dds), name = name.formula[i]))
                    }, 
                    mc.cores = 10)
    
    eFDR <- mclapply(1:length(vars.pool),
                     function(i) qvalue(res[[i]][["pvalue"]])$qvalues,
                     mc.cores = 10)
    
    eFDR.GeneName <- mclapply(1:length(vars.pool),
                              function(i) {
                                q.val <- eFDR[[i]]
                                q.ord <- order(q.val)[1:threshold.val]
                                q.val <- q.val[q.ord]
                                names(q.val) <- rownames(cts)[q.ord]
                                return(q.val)
                              },
                              mc.cores = 10)
    
    names(eFDR.GeneName) = vars.pool
    
    return(eFDR.GeneName)
  }
  
  eFDR.GeneName_train <- eFDR.GeneName_test <- NULL
  
  if (use.train) {
    result_train <- run_DESeq2(as.matrix(dat.expr.train), dat.pheno.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
  } else {
    result_train <- run_DESeq2(as.matrix(dat.expr.train), dat.pheno.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
    result_test <- run_DESeq2(as.matrix(dat.expr.test), dat.pheno.test, 5000)
    eFDR.GeneName_test <- result_test$eFDR.GeneName
  }
}

multi_TMM_top_edgeR <- function(seed.val = 999, use.train = TRUE) {
  run_edgeR <- function(dat.expr, design, threshold.val) {
    y <- DGEList(counts = dat.expr, lib.size = colSums(dat.expr)) |>
      calcNormFactors() |>
      estimateGLMTrendedDisp(design)
    
    edgeR.fit <- glmQLFit(y, design)
    
    qlf <- mclapply(1:length(vars.pool),
                    function(i) glmQLFTest(edgeR.fit, coef = i + 1),
                    mc.cores = 10)
    
    eFDR <- mclapply(1:length(vars.pool),
                     function(i) qlf[[i]][["table"]][["PValue"]] |>
                       qvalue(),
                     mc.cores = 10)
    
    eFDR.GeneName <- mclapply(1:length(vars.pool),
                              function(i) {
                                q.val <- eFDR[[i]][["qvalues"]]
                                q.ord <- order(q.val)[1:threshold.val]
                                q.val <- q.val[q.ord]
                                names(q.val) <- rownames(dat.expr)[q.ord]
                                return(q.val)
                              },
                              mc.cores = 10)
    
    names(eFDR.GeneName) = vars.pool
    
    return(eFDR.GeneName)
  }
  
  eFDR.GeneName_train <- eFDR.GeneName_test <- NULL
  
  if (use.train) {
    result_train <- run_edgeR(dat.expr.train, design.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
  } else {
    result_train <- run_edgeR(dat.expr.train, design.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
    result_test <- run_edgeR(dat.expr.test, design.test, 5000)
    eFDR.GeneName_test <- result_test$eFDR.GeneName
  }
}

multi_TMM_top_eBayes <- function(seed.val = 999, use.train = TRUE) {
  run_eBayes <- function(dat.expr, design, threshold.val) {
    voom.data <- voom(counts = dat.expr,
                      design = design,
                      lib.size = colSums(dat.expr) * calcNormFactors(dat.expr))
    
    eBayes.fit <- lmFit(voom.data[["E"]],
                        design = design,
                        weights = voom.data[["weights"]]) |> eBayes()
    
    eFDR <- mclapply(1:length(vars.pool),
                     function(i) eBayes.fit[["p.value"]][, i + 1] |> qvalue(),
                     mc.cores = 10)
    
    eFDR.GeneName <- mclapply(1:length(vars.pool),
                              function(i) {
                                q.val <- eFDR[[i]][["qvalues"]]
                                return(q.val[order(q.val)[1:threshold.val]])
                              },
                              mc.cores = 10)
    
    names(eFDR.GeneName) = vars.pool
    
    return(eFDR.GeneName)
  }
  
  eFDR.GeneName_train <- eFDR.GeneName_test <- NULL
  
  if (use.train) {
    result_train <- run_eBayes(dat.expr.train, design.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
  } else {
    result_train <- run_eBayes(dat.expr.train, design.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
    result_test <- run_eBayes(dat.expr.test, design.test, 5000)
    eFDR.GeneName_test <- result_test$eFDR.GeneName
  }
}

multi_TMM_top_voom.limma <- function(seed.val = 999, use.train = TRUE) {
  run_voom <- function(dat.expr, design, threshold.val) {
    voom.data <- voom(counts = dat.expr,
                      design = design,
                      lib.size = colSums(dat.expr) * calcNormFactors(dat.expr))
    
    voom.fit <- lmFit(voom.data[["E"]],
                      design = design,
                      weights = voom.data[["weights"]])
    
    eFDR <- mclapply(1:length(vars.pool),
                     function(i) {
                       t <- voom.fit[["coefficients"]][, i + 1] / voom.fit[["stdev.unscaled"]][, i + 1] / voom.fit[["sigma"]]
                       p <- 2 * pt(-abs(t), df = voom.fit[["df.residual"]])
                       return(qvalue(p))
                     },
                     mc.cores = 10)
    
    eFDR.GeneName <- mclapply(1:length(vars.pool),
                              function(i) {
                                q.val <- eFDR[[i]][["qvalues"]]
                                return(q.val[order(q.val)[1:threshold.val]])
                              },
                              mc.cores = 10)
    
    names(eFDR.GeneName) = vars.pool
    
    return(eFDR.GeneName)
  }
  
  eFDR.GeneName_train <- eFDR.GeneName_test <- NULL
  
  if (use.train) {
    result_train <- run_voom(dat.expr.train, design.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
  } else {
    result_train <- run_voom(dat.expr.train, design.train, 5000)
    eFDR.GeneName_train <- result_train$eFDR.GeneName
    result_test <- run_voom(dat.expr.test, design.test, 5000)
    eFDR.GeneName_test <- result_test$eFDR.GeneName
  }
}
