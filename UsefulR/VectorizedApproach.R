t0 <- microbenchmark::get_nanotime()
no.cDEGs <- length(cDEGs.all)
cDEGs.mat <- matrix(0, nrow = no.cDEGs, ncol = 10) 
rownames(cDEGs.mat) <- cDEGs.all
colnames(cDEGs.mat) <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM",
                         "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM", "Class")

# Conventional for loop
for (i in 1:no.cDEGs) { 
  if (cDEGs.all[i] %in% BMAseq.cDEGs) { 
    cDEGs.mat[i, 1] <- 1
  } 
  if (cDEGs.all[i] %in% DESeq2.UVM.cDEGs) {
    cDEGs.mat[i, 2] <- 1
  } 
  if (cDEGs.all[i] %in% DESeq2.MVM.cDEGs) {
    cDEGs.mat[i, 3] <- 1
  } 
  if (cDEGs.all[i] %in% edgeR.UVM.cDEGs) {
    cDEGs.mat[i, 4] <- 1
  } 
  if (cDEGs.all[i] %in% edgeR.MVM.cDEGs) {
    cDEGs.mat[i, 5] <- 1
  } 
  if (cDEGs.all[i] %in% eBayes.UVM.cDEGs) {
    cDEGs.mat[i, 6] <- 1
  } 
  if (cDEGs.all[i] %in% eBayes.MVM.cDEGs) {
    cDEGs.mat[i, 7] <- 1
  } 
  if (cDEGs.all[i] %in% voom.limma.UVM.cDEGs) {
    cDEGs.mat[i, 8] <- 1
  } 
  if (cDEGs.all[i] %in% voom.limma.MVM.cDEGs) {
    cDEGs.mat[i, 9] <- 1
  }
  cDEGs.mat[i, 10] <- paste(cDEGs.mat[i, 1:9], collapse = "")
}
t1 <- microbenchmark::get_nanotime()
t1 - t0

# 2x faster using vectorization
t0 <- microbenchmark::get_nanotime()
no.cDEGs <- length(cDEGs.all)
cDEGs.mat.2 <- matrix(0, nrow = no.cDEGs, ncol = 10) 
rownames(cDEGs.mat.2) <- cDEGs.all
colnames(cDEGs.mat.2) <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM",
                           "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM", "Class")
cDEGs.mat.2[, 1] <- as.integer(cDEGs.all %in% BMAseq.cDEGs)
cDEGs.mat.2[, 2] <- as.integer(cDEGs.all %in% DESeq2.UVM.cDEGs)
cDEGs.mat.2[, 3] <- as.integer(cDEGs.all %in% DESeq2.MVM.cDEGs)
cDEGs.mat.2[, 4] <- as.integer(cDEGs.all %in% edgeR.UVM.cDEGs)
cDEGs.mat.2[, 5] <- as.integer(cDEGs.all %in% edgeR.MVM.cDEGs)
cDEGs.mat.2[, 6] <- as.integer(cDEGs.all %in% eBayes.UVM.cDEGs)
cDEGs.mat.2[, 7] <- as.integer(cDEGs.all %in% eBayes.MVM.cDEGs)
cDEGs.mat.2[, 8] <- as.integer(cDEGs.all %in% voom.limma.UVM.cDEGs)
cDEGs.mat.2[, 9] <- as.integer(cDEGs.all %in% voom.limma.MVM.cDEGs)
cDEGs.mat.2[, 10] <- apply(cDEGs.mat.2[, 1:9], 1, function(x) paste(x, collapse = ""))
t1 <- microbenchmark::get_nanotime()
t1 - t0
