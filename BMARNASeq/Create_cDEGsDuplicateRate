cDEGs <- NULL
threshold.vec <- seq(1000, 5000, 1000)
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)

duplicated_rate_cDEGs <- function(var.name = "BMI", method.name = "BMAseq") {
  for (threshold.i in threshold.vec) {
    for (seed.i in seed.vec) {
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sMultiInt%s.RData", method.name, seed.i))
      cDEGs <- c(cDEGs, intersect(names(get(sprintf("%s.eFDR.Main.train", method.name))[[var.name]][1:threshold.i]), 
                                  names(get(sprintf("%s.eFDR.Main.test", method.name))[[var.name]][1:threshold.i])))
    }
    print(paste("The duplicated rate of cDEGs associated with", var.name, "identified by", method.name, "among 10 random seeds for the threshold", threshold.i, "is", sum(duplicated(cDEGs))/length(cDEGs)))
    cDEGs <- NULL
  }
}

mapply(duplicated_rate_cDEGs,
       c("BMI","AGE","SEX","WBC","BMIxSEX"), 
       c("BMAseq", "DESeq2", "edgeR", "eBayes", "VoomLimma"))
