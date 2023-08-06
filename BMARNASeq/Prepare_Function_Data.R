##------Attach libraries and functions------
# devtools::install("/Library/Frameworks/R.framework/Versions/4.2/Resources/library/BMAseq") {BMAseq} (2023Apr1 version)  
easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", 
                        "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", 
                        "openxlsx", "grDevices", "viridis", "cowplot", "clusterProfiler",
                        "org.Hs.eg.db", "latex2exp", "enrichplot") |> suppressPackageStartupMessages()

Bayesfactor <- BMAseq:::Bayesfactor

##------Process data------
dat.expr <- dget("../ApplicationData/derived/dat.expr.Subcutaneous") 
dat.pheno <- dget("../ApplicationData/derived/dat.pheno.Subcutaneous") 
dim(dat.expr) # 24660 genes and 404 subjects
dim(dat.pheno) # 404 subjects and 13 phenotypes
dat.pheno[1:5, ]
dat.expr[1:5, 1:5]
paste0("The column names of gene expression data", ifelse(all(colnames(dat.expr) == rownames(dat.pheno)), " MATCH ", "NOT MATCH"), "the row names of phenotype data.")

##------Pre-filter the genes------
# Here we perform the median absolute deviation with the threshold of 0.8 to select genes that are most likely to distinguish the samples
dat.expr.new <- dat.expr[apply(dat.expr, 1, function(x) mad(x) > 0.8), ] # We have 24455 genes

seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
for (seed.i in seed.vec) {
  
  # Divide the gene expression data and phenotype data into the 50% training set and 50% test set 
  test.ind <- sample(1:nrow(dat.pheno), ceiling(0.5*nrow(dat.pheno)))
  dat.pheno.train <- dat.pheno[-test.ind, ]
  dat.pheno.test <- dat.pheno[test.ind, ]
  dat.expr.train <- dat.expr.new[, rownames(dat.pheno.train)]
  dat.expr.test <- dat.expr.new[, rownames(dat.pheno.test)]
  
  # Check if the gene expression data matches the phenotype data  
  print(paste0("The column names of training gene expression data", ifelse(all(colnames(dat.expr.train) == rownames(dat.pheno.train)), " MATCH ", "NOT MATCH"), "the row names of training phenotype data."))
  print(paste0("The column names of test gene expression data", ifelse(all(colnames(dat.expr.test) == rownames(dat.pheno.test)), " MATCH ", "NOT MATCH"), "the row names of test phenotype data."))
  
  # Make the table one of patient characteristics per random seed trial
  dat.pheno.train$group <- "Training"
  dat.pheno.test$group <- "Test"
  dat.pheno.all <- rbind(dat.pheno.train, dat.pheno.test)
  tab1 <- tbl_summary(data = dat.pheno.all,
                      by = "group") %>% 
    add_p() %>% 
    as_gt() %>%
    gt::gtsave(filename = paste0("../ApplicationResult/Multi/RandomSeed/TMM_Top2000/", date.analysis, "_", "Table1", "_", seed.num, ".rtf"))
  
  # Save gene expression and phenotype data
  saveRDS(dat.expr.train, 
          file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i))
  saveRDS(dat.expr.test, 
          file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i))
  saveRDS(dat.pheno.train, 
          file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i))
  saveRDS(dat.pheno.test, 
          file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i))
}


##------Modify `BMAseq.multi.postprob` to be compatible with the trimmed-mean of M-values (TMM)------
BMAseq.multi.postprob.norm <- function(dat.expr.counts, dat.pheno, var.pool, max.nvar, interaction = NULL, cut.BF = 1) {
  
  
  if (sum(colnames(dat.expr.counts) != rownames(dat.pheno)) > 0)
    stop("Two datasets dat.expr.counts and dat.pheno must be matched by subjects.")
  if (is.null(var.pool) || length(var.pool) < 2)
    stop("Must provide at least two variables of interest.")
  if (sum(var.pool %in% colnames(dat.pheno)) != length(var.pool))
    stop("Variables of interest must be in datasets dat.pheno.")
  if (is.null(max.nvar))
    stop("Must provide max.nvar.")
  
  
  n.var <- length(var.pool)
  n.gene <- nrow(dat.expr.counts)
  
  
  # estimate voom weights
  y0 <- voom(counts = dat.expr.counts, 
             lib.size = colSums(dat.expr.counts)*calcNormFactors(dat.expr.counts))
  input.dat.expr <- y0$E
  input.weights <- y0$weights
  
  
  # build model space
  result.modelspace <- Modelspace(dat.pheno, var.pool, max.nvar, interaction = interaction)
  input.model.space <- result.modelspace$model.space
  input.dat.pheno <- result.modelspace$dat.pheno.new
  
  
  # calculate Bayes factor
  out.bf <- vapply(1:length(input.model.space), function(i) {
    print(paste0(i, ". ", input.model.space[i]))
    design <- model.matrix(formula(input.model.space[i]), data = input.dat.pheno)
    colnames(design)[1] <- "Int"
    vapply(1:n.gene, function(k) Bayesfactor(design, input.dat.expr[k, ], input.weights[k, ]), FUN.VALUE = as.double(1))
  }, FUN.VALUE = as.double(1:n.gene))
  
  
  # obtain prior model probability
  bf.max <- t(apply(out.bf, 1, function(x) ifelse(x == max(x) & x > cut.BF, 1, 0)))
  prior.modelprob <- apply(bf.max, 2, sum)/n.gene
  pm.est <- prior.modelprob
  n.model <- length(prior.modelprob)
  pm.prior0 <- pm.prior1 <- rep(0.5/(n.model - 1), n.model - 1)
  alpha <- 0.2
  for (j in 1:30) {
    pm.prior1 <- pm.est[-1]/apply(t(apply(out.bf[, -1], 1, function(x) x/(1 - sum(pm.prior0) + sum(x * pm.prior0)))), 2,
                                  mean)
    pm.prior0 <- pm.prior0 + (pm.prior1 - pm.prior0) * alpha
  }
  prior.modelprob <- c(1 - sum(pm.prior0), pm.prior0)
  
  
  # calculate posterior model probability
  post.modelprob <- t(apply(out.bf, 1, function(x) x * prior.modelprob/sum(x * prior.modelprob)))
  rownames(post.modelprob) <- rownames(dat.expr.counts)
  colnames(post.modelprob) <- input.model.space
  
  
  
  return(list(dat.expr.logcpm = input.dat.expr, weights = input.weights, model.space = input.model.space, dat.pheno.new = input.dat.pheno,
              post.modelprob = post.modelprob, var.pool = var.pool, interaction = interaction))
}
