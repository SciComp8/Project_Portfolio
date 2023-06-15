model.type.vec <- c("Multi", "Uni")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
model.method.vec <- c("BMAseq", "DESeq2", "edgeR", "eBayes", "VoomLimma")

# Test
# seed.i <- 8809678
# model.method.i <- "VoomLimma"
# model.type.i <- "Uni"

for (model.type.i in model.type.vec) {
  for (seed.i in seed.vec) {
    for (model.method.i in model.method.vec) {
      if (model.type.i == "Multi" & model.method.i == "BMAseq") {
        file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/%s%s%s.RData", 
                             model.type.i, model.method.i, model.type.i, seed.i)
        data.env <- local({load(file.name); environment()})
        names(data.env$BMAseq.eFDR.Main.train) = names(data.env$BMAseq.eFDR.Main.test) = var.vec
        save(BMAseq.eFDR.Main.train, BMAseq.eFDR.Main.test, 
             envir = data.env,
             file = sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                            model.type.i, model.method.i, model.type.i, seed.i))
        rm(list = "data.env")
      } else if (model.type.i == "Uni" & model.method.i == "BMAseq") {
        next # skip the current iteration and start next iteration
      } else if (model.method.i == "DESeq2") {
        file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/%s%s%s.RData", 
                             model.type.i, model.method.i, model.type.i, seed.i)
        data.env <- local({load(file.name); environment()})
        save(DESeq2.eFDR.GeneName.train, DESeq2.eFDR.GeneName.test, 
             envir = data.env,
             file = sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                            model.type.i, model.method.i, model.type.i, seed.i))
        rm(list = "data.env")
      } else if (model.method.i == "edgeR") {
        file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/%s%s%s.RData", 
                             model.type.i, model.method.i, model.type.i, seed.i)
        data.env <- local({load(file.name); environment()})
        save(edgeR.eFDR.GeneName.train, edgeR.eFDR.GeneName.test, 
             envir = data.env,
             file = sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                            model.type.i, model.method.i, model.type.i, seed.i))
        rm(list = "data.env")
      } else if (model.method.i == "eBayes") {
        file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/%s%s%s.RData", 
                             model.type.i, model.method.i, model.type.i, seed.i)
        data.env <- local({load(file.name); environment()})
        save(eBayes.eFDR.train2, eBayes.eFDR.test2, 
             envir = data.env,
             file = sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                            model.type.i, model.method.i, model.type.i, seed.i))
        rm(list = "data.env")
      } else if (model.method.i == "VoomLimma") {
        file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/%s%s%s.RData", 
                             model.type.i, model.method.i, model.type.i, seed.i)
        data.env <- local({load(file.name); environment()})
        save(voom.eFDR.train2, voom.eFDR.test2, 
             envir = data.env,
             file = sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                            model.type.i, model.method.i, model.type.i, seed.i))
        rm(list = "data.env")
      }
    }
  }
}
