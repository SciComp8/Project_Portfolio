# ! Eliminate code repetition by encapsulating the common DESeq2 operations

# Redundant code with 112 lines
multi_TMM_top_DESeq2_old <- function(seed.val = 999, threshold = 2000, use_train = TRUE) {
  cts.train <- as.matrix(dat.expr.train)
  cts.test <- as.matrix(dat.expr.test)
  coldata.train <- dat.pheno.train
  coldata.test <- dat.pheno.test
  
  name.formula <- c("BMI_high_vs_low", "AGE_old_vs_young", "SEX_male_vs_female", "MHABNWBC_yes_vs_no") 
  
  if (use_train) {
    lib.size <- colSums(cts.train)
    norm.factor <- calcNormFactors(cts.train, method = "TMM")
    size.factor <- lib.size*norm.factor/exp(mean(log(lib.size*norm.factor))) 
    
    dds <- DESeqDataSetFromMatrix(countData = cts.train, 
                                  colData = coldata.train, 
                                  design = ~BMI + AGE + SEX + MHABNWBC)
    sizeFactors(dds) <- size.factor
    
    res.train <- mclapply(1:length(vars.pool), 
                          function(i) {
                            return(results(DESeq(dds), name = name.formula[i])) },
                          mc.cores = 10)
    
    DESeq2.eFDR.train <- mclapply(1:length(vars.pool),
                                  function(i) qvalue(res.train[[i]][["pvalue"]])$qvalues,
                                  mc.cores = 10)
    
    DESeq2.eFDR.train2 <- mclapply(1:length(vars.pool),
                                   function(i) {
                                     q.val <- DESeq2.eFDR.train[[i]]
                                     return(q.val[order(q.val)[1:threshold]])
                                   },
                                   mc.cores = 10)
    
    DESeq2.eFDR.GeneName.train <- mclapply(1:length(vars.pool),
                                           function(i) {
                                             q.val <- DESeq2.eFDR.train[[i]]
                                             return(rownames(cts.train)[order(q.val)[1:threshold]])
                                           },
                                           mc.cores = 10)
    
    names(res.train) = names(DESeq2.eFDR.train) = names(DESeq2.eFDR.train2) = names(DESeq2.eFDR.GeneName.train) = vars.pool 
  } else {
    lib.size <- colSums(cts.train)
    norm.factor <- calcNormFactors(cts.train, method = "TMM")
    size.factor <- lib.size*norm.factor/exp(mean(log(lib.size*norm.factor))) 
    
    dds <- DESeqDataSetFromMatrix(countData = cts.train, 
                                  colData = coldata.train, 
                                  design = ~BMI + AGE + SEX + MHABNWBC)
    sizeFactors(dds) <- size.factor
    
    res.train <- mclapply(1:length(vars.pool), 
                          function(i) {
                            return(results(DESeq(dds), name = name.formula[i])) },
                          mc.cores = 10)
    
    DESeq2.eFDR.train <- mclapply(1:length(vars.pool),
                                  function(i) qvalue(res.train[[i]][["pvalue"]])$qvalues,
                                  mc.cores = 10)
    
    DESeq2.eFDR.train2 <- mclapply(1:length(vars.pool),
                                   function(i) {
                                     q.val <- DESeq2.eFDR.train[[i]]
                                     return(q.val[order(q.val)[1:threshold]])
                                   },
                                   mc.cores = 10)
    
    DESeq2.eFDR.GeneName.train <- mclapply(1:length(vars.pool),
                                           function(i) {
                                             q.val <- DESeq2.eFDR.train[[i]]
                                             return(rownames(cts.train)[order(q.val)[1:threshold]])
                                           },
                                           mc.cores = 10)
    
    names(res.train) = names(DESeq2.eFDR.train) = names(DESeq2.eFDR.train2) = names(DESeq2.eFDR.GeneName.train) = vars.pool 
    
    lib.size <- colSums(cts.test)
    norm.factor <- calcNormFactors(cts.test, method = "TMM")
    size.factor <- lib.size*norm.factor/exp(mean(log(lib.size*norm.factor))) 
    
    dds <- DESeqDataSetFromMatrix(countData = cts.test, 
                                  colData = coldata.test, 
                                  design = ~BMI + AGE + SEX + MHABNWBC)
    sizeFactors(dds) <- size.factor
    
    res.test <- mclapply(1:length(vars.pool), 
                         function(i) {
                           return(results(DESeq(dds), name = name.formula[i])) }, 
                         mc.cores = 10)
    
    DESeq2.eFDR.test <- mclapply(1:length(vars.pool),
                                 function(i) qvalue(res.test[[i]][["pvalue"]])$qvalues,
                                 mc.cores = 10)
    
    DESeq2.eFDR.test2 <- mclapply(1:length(vars.pool),
                                  function(i) {
                                    q.val <- DESeq2.eFDR.test[[i]]
                                    return(q.val[order(q.val)[1:threshold]])
                                  },
                                  mc.cores = 10)
    
    DESeq2.eFDR.GeneName.test <- mclapply(1:length(vars.pool),
                                          function(i) {
                                            q.val <- DESeq2.eFDR.test[[i]]
                                            return(rownames(cts.test)[order(q.val)[1:threshold]])
                                          },
                                          mc.cores = 10)
    
    names(res.test) = names(DESeq2.eFDR.test) = names(DESeq2.eFDR.test2) = names(DESeq2.eFDR.GeneName.test) = vars.pool
  }
}

# Concise code with 63 lines
run_DESeq2 <- function(cts, coldata, threshold) {
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
  
  eFDR2 <- mclapply(1:length(vars.pool),
                    function(i) {
                      q.val <- eFDR[[i]]
                      return(q.val[order(q.val)[1:threshold]])
                    },
                    mc.cores = 10)
  
  GeneName <- mclapply(1:length(vars.pool),
                       function(i) {
                         q.val <- eFDR[[i]]
                         return(rownames(cts)[order(q.val)[1:threshold]])
                       },
                       mc.cores = 10)
  
  names(res) = names(eFDR) = names(eFDR2) = names(GeneName) = vars.pool
  
  return(list(res = res, eFDR = eFDR, eFDR2 = eFDR2, GeneName = GeneName))
}

multi_TMM_top_DESeq2 <- function(seed.val = 999, threshold = 2000, use_train = TRUE) {
  results_train <- eFDR_train <- eFDR2_train <- GeneName_train <- NULL
  results_test <- eFDR_test <- eFDR2_test <- GeneName_test <- NULL
  
  if (use_train) {
    result_train <- run_DESeq2(as.matrix(dat.expr.train), dat.pheno.train, threshold)
    results_train <- result_train$res
    eFDR_train <- result_train$eFDR
    eFDR2_train <- result_train$eFDR2
    GeneName_train <- result_train$GeneName
  } else {
    result_train <- run_DESeq2(as.matrix(dat.expr.train), dat.pheno.train, threshold)
    results_train <- result_train$res
    eFDR_train <- result_train$eFDR
    eFDR2_train <- result_train$eFDR2
    GeneName_train <- result_train$GeneName
    
    result_test <- run_DESeq2(as.matrix(dat.expr.test), dat.pheno.test, threshold)
    results_test <- result_test$res
    eFDR_test <- result_test$eFDR
    eFDR2_test <- result_test$eFDR2
    GeneName_test <- result_test$GeneName
  }
}


