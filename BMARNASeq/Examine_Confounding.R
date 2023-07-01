easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", "openxlsx", "cowplot") |> 
  suppressPackageStartupMessages()
Bayesfactor <- BMAseq:::Bayesfactor

# Load data
load("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti8809678.RData")

# BMI ~ AGE
with(dat.pheno.train, chisq.test(BMI, AGE))$p.value
with(dat.pheno.test, chisq.test(BMI, AGE))$p.value

# BMI ~ each gene expression
BMI.train <- dat.pheno.train$BMI |> factor() |> as.numeric()
BMI.test <- dat.pheno.test$BMI |> factor() |> as.numeric()
p.train <-
  sapply(1:nrow(dat.expr.train), function(x) {
    wilcox.test(dat.expr.train[x, ] |> c() |> unlist(),
                BMI.train)$p.value
  })
p.test <-
  sapply(1:nrow(dat.expr.test), function(x) {
    wilcox.test(dat.expr.test[x, ] |> c() |> unlist(),
                BMI.test)$p.value
  })
BMI.p.df <- data.frame(gene.id = rownames(dat.expr.train), p.train = p.train, p.test = p.test)


# AGE ~ each gene expression
AGE.train <- dat.pheno.train$AGE |> factor() |> as.numeric()
AGE.test <- dat.pheno.test$AGE |> factor() |> as.numeric()
p.train <-
  sapply(1:nrow(dat.expr.train), function(x) {
    wilcox.test(dat.expr.train[x, ] |> c() |> unlist(),
                AGE.train)$p.value
  })
p.test <-
  sapply(1:nrow(dat.expr.test), function(x) {
    wilcox.test(dat.expr.test[x, ] |> c() |> unlist(),
                AGE.test)$p.value
  })
AGE.p.df <- data.frame(gene.id = rownames(dat.expr.train), p.train = p.train, p.test = p.test)
