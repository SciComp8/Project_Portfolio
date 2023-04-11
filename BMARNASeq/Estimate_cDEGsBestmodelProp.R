library(data.table)
load("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti120.RData") 
View(output.multi.train.s2)
temp <- rbind(output.multi.train.s2[["DEG.bestmodel.Main"]][["BMI"]], output.multi.test.s2[["DEG.bestmodel.Main"]][["BMI"]])
cDEG.all <- intersect(BMAseq.eFDR.Main.train$BMI[1:2000] |> names(), BMAseq.eFDR.Main.test$BMI[1:2000] |> names())  # ranking threshold: 2000
cDEG.bestmodel <- temp[temp[, 2] %in% cDEG.all, ]
num.cDEG.bestmodel <- nrow(cDEG.bestmodel)
prop.bestunimodel <- sum(cDEG.bestmodel[, 3] == "~1+BMI") / num.cDEG.bestmodel
prop.bestmultimodel <- 1 - prop.bestunimodel - sum(cDEG.bestmodel[, 3] == "~1") / num.cDEG.bestmodel
prop.dt <- data.table(threshold = 2000, seed = 120, variable = "BMI", num.cDEG.bestmodel = num.cDEG.bestmodel, prop.bestunimodel = prop.bestunimodel, prop.bestmultimodel = prop.bestmultimodel)

calc_prop_bestmodel <- function (threshold = NULL, seed = NULL, var.name = NULL) {
  temp <- rbind(output.multi.train.s2[["DEG.bestmodel.Main"]][[var.name]], output.multi.test.s2[["DEG.bestmodel.Main"]][[var.name]])
  cDEG.all <- intersect(BMAseq.eFDR.Main.train[[var.name]][1:threshold] |> names(), BMAseq.eFDR.Main.test[[var.name]][1:threshold] |> names())
  cDEG.bestmodel <- temp[temp[, 2] %in% cDEG.all, ]
  num.cDEG.bestmodel <- nrow(cDEG.bestmodel)
  prop.bestunimodel <- sum(cDEG.bestmodel[, 3] == sprintf("~1+%s", var.name)) / num.cDEG.bestmodel
  prop.bestmultimodel <- 1 - prop.bestunimodel - sum(cDEG.bestmodel[, 3] == "~1") / num.cDEG.bestmodel
  prop.dt <- data.table(threshold = threshold, seed = seed, variable = var.name, num.cDEG.bestmodel = num.cDEG.bestmodel, prop.bestunimodel = prop.bestunimodel, prop.bestmultimodel = prop.bestmultimodel)
  return(prop.dt)
}

prop.list <- NULL
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- seq(1000, 5000, 1000)
for (seed.i in seed.vec) {
  load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%s.RData", seed.i))
  for (threshold.i in threshold.vec) {
    for (var.i in var.vec) {
      prop.i <- calc_prop_bestmodel(threshold = threshold.i, seed = seed.i, var.name = var.i)
      prop.list <- append(prop.list, list(prop.i))
    }
  }
}

saveRDS(rbindlist(prop.list), file = "../ApplicationResult/Multi/RandomSeed/BestModelProportion/BMAseq.bestmodel.prop.RDS")
BMAseq.bestmodel.prop <- readRDS("../ApplicationResult/Multi/RandomSeed/BestModelProportion/BMAseq.bestmodel.prop.RDS")

