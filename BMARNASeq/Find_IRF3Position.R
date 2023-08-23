##------IRF3 position in the ranking list------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
best.model <- best.model |> column_to_rownames(var = "...1")
# SYMBOL         ENSEMBL
# 1   IRF3 ENSG00000126456

# Find which random seed trial shows IRF3 as the cDEG detected with BMAseq at a 5000 ranking threshold
IRF3.best.model <- best.model[rownames(best.model) == "ENSG00000126456.11", 1:20]
IRF3.best.model <- IRF3.best.model |> unlist() # Transform the dataframe into a named vector
seed.trial <- character()
for (name in names(IRF3.best.model)) {
  if (!is.na(IRF3.best.model[name])) {
    seed.trial.i <- sub("bestmodel.", "", name)
    seed.trial <- c(seed.trial, seed.trial.i)
  }
}

var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
IRF3.pos <- numeric()
for (seed.trial.i in seed.trial) {
  seed <- sub(".*\\.", "", seed.trial.i)
  data.type <- sub("\\..*", "", seed.trial.i)
  temp <- new.env()
  target.file <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%s.RData", seed)
  load(target.file, envir = temp)
  names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
  IRF3.pos.i <- which(names(temp[[sprintf("BMAseq.eFDR.Main.%s", data.type)]][["BMI"]]) == "ENSG00000126456.11")
  IRF3.pos <- c(IRF3.pos, IRF3.pos.i)
  rm("temp")
}

names(IRF3.pos) <- seed.trial
IRF3.pos |> data.frame()
