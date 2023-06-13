easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", "openxlsx") |> 
  suppressPackageStartupMessages()

threshold <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")

get_venndata <- function(seed.i = 8809678, 
                         var.name = "BMI", 
                         model.type.1 = "Multi", 
                         model.type.2 = "Uni",
                         model.method.1 = "BMAseq",
                         model.method.2 = "DESeq2",
                         threshold = 5000,
                         ...){
  file.name.1 <- sprintf("../ApplicationData/derived/RandomSeed/Top%s/%sModel/%s%s%s.RData", 
                         threshold, model.type.1, model.method.1, model.type.1, seed.i)
  load(file.name.1)
  if (model.method.1 == "BMAseq") {
    names(BMAseq.eFDR.Main.train) = names(BMAseq.eFDR.Main.test) = var.vec
    cDEGs.1 <- intersect(names(BMAseq.eFDR.Main.train[[var.name]][1:threshold]), 
                         names(BMAseq.eFDR.Main.test[[var.name]][1:threshold]))
  } else if (model.method.1 == "DESeq2") {
    cDEGs.1 <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                         DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.1 == "edgeR") {
    cDEGs.1 <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                         edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.1 == "eBayes") {
    cDEGs.1 <- intersect(names(edgeR.eFDR.GeneName.train[[var.name]][1:threshold]),
                         names(edgeR.eFDR.GeneName.test[[var.name]][1:threshold]))
  } else if (model.method.1 == "VoomLimma") {
    cDEGs.1 <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold]),
                         names(voom.eFDR.train2[[var.name]][1:threshold]))
  }
  file.name.2 <- sprintf("../ApplicationData/derived/RandomSeed/Top%s/%sModel/%s%s%s.RData", 
                         threshold, model.type.2, model.method.2, model.type.2, seed.i)
  load(file.name.2)
  if (model.method.2 == "DESeq2") {
    cDEGs.2 <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                         DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.2 == "edgeR") {
    cDEGs.2 <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                         edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.2 == "eBayes") {
    cDEGs.2 <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold]),
                         names(eBayes.eFDR.train2[[var.name]][1:threshold]))
  } else if (model.method.2 == "VoomLimma") {
    cDEGs.2 <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold]),
                         names(voom.eFDR.train2[[var.name]][1:threshold]))
  }
  return(list(cDEGs.1, cDEGs.2))
}
