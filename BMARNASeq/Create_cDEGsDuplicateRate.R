cDEGs <- NULL
threshold.vec <- seq(1000, 5000, 1000)
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)

duplicated_rate_cDEGs <- function(var.name = "BMI", method.name = "BMAseq") {
  for (threshold.i in threshold.vec) {
    for (seed.i in seed.vec) {
      if (method.name == "BMAseq" & var.name != "BMIxSEX") {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sMultiInt%s.RData", method.name, seed.i))
        cDEGs <- c(cDEGs, intersect(names(get(sprintf("%s.eFDR.Main.train", method.name))[[var.name]][1:threshold.i]), 
                                    names(get(sprintf("%s.eFDR.Main.test", method.name))[[var.name]][1:threshold.i])))
        } else if (method.name == "BMAseq" & var.name == "BMIxSEX") {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sMultiInt%s.RData", method.name, seed.i))
        var.name2 <- "BMI"
        cDEGs <- c(cDEGs, intersect(names(get(sprintf("%s.eFDR.Interaction.train", method.name))[[var.name2]][1:threshold.i]), 
                                    names(get(sprintf("%s.eFDR.Interaction.test", method.name))[[var.name2]][1:threshold.i])))
      } else if (method.name %in% c("DESeq2", "edgeR")) {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sMultiInt%s.RData", method.name, seed.i))
        cDEGs <- c(cDEGs, intersect(get(sprintf("%s.eFDR.GeneName.train", method.name))[[var.name]][1:threshold.i], 
                                    get(sprintf("%s.eFDR.GeneName.test", method.name))[[var.name]][1:threshold.i]))
      } else if (method.name == "VoomLimma") {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sMultiInt%s.RData", method.name, seed.i))
        method.name2 <- "voom"
        cDEGs <- c(cDEGs, intersect(names(get(sprintf("%s.eFDR.train2", method.name2))[[var.name]][1:threshold.i]), 
                                    names(get(sprintf("%s.eFDR.test2", method.name2))[[var.name]][1:threshold.i])))
      } else {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sMultiInt%s.RData", method.name, seed.i))
        cDEGs <- c(cDEGs, intersect(names(get(sprintf("%s.eFDR.train2", method.name))[[var.name]][1:threshold.i]), 
                                    names(get(sprintf("%s.eFDR.test2", method.name))[[var.name]][1:threshold.i])))
      }
    }
    print(paste("The duplicated rate of cDEGs associated with", var.name, "identified by", method.name, "among 10 random seeds for the threshold", threshold.i, "is", sum(duplicated(cDEGs))/length(cDEGs)))
    cDEGs <- NULL
  }
}
