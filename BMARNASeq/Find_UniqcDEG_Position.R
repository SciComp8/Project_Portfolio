##------IRF3 position in the ranking list------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
best.model <- best.model |> column_to_rownames(var = "...1")
# SYMBOL         ENSEMBL
# 1   IRF3 ENSG00000126456

# Special case for one seed
temp <- new.env()
target.file <- "../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti8809678.RData"
load(target.file, envir = temp)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
IRF3.pos.train <- which(names(temp$BMAseq.eFDR.Main.train[["BMI"]]) == "ENSG00000126456.11")
IRF3.pos.test <- which(names(temp$BMAseq.eFDR.Main.test[["BMI"]]) == "ENSG00000126456.11")

# General case for all possible seeds
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
IRF3.pos |> data.frame() |> View()


##------DDX1 position in the ranking list------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
best.model <- best.model |> column_to_rownames(var = "...1")
# SYMBOL         ENSEMBL
# 1   DDX1 ENSG00000079785

# Special case for one seed
temp <- new.env()
target.file <- "../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti8809678.RData"
load(target.file, envir = temp)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
DDX1.pos.train <- which(grepl(pattern = "ENSG00000079785", x = names(temp$BMAseq.eFDR.Main.train[["BMI"]])))
DDX1.pos.test <- which(grepl(pattern = "ENSG00000079785", x = names(temp$BMAseq.eFDR.Main.test[["BMI"]])))

# General case for all possible seeds
# Find which random seed trial shows DDX1 as the cDEG detected with BMAseq at a 5000 ranking threshold
DDX1.best.model <- best.model[which(grepl(pattern = "ENSG00000079785", x = rownames(best.model))), 1:20]
DDX1.best.model <- DDX1.best.model |> unlist() # Transform the dataframe into a named vector
seed.trial <- character()
for (name in names(DDX1.best.model)) {
  if (!is.na(DDX1.best.model[name])) {
    seed.trial.i <- sub("bestmodel.", "", name)
    seed.trial <- c(seed.trial, seed.trial.i)
  }
}

var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
DDX1.pos <- numeric()
for (seed.trial.i in seed.trial) {
  seed <- sub(".*\\.", "", seed.trial.i)
  data.type <- sub("\\..*", "", seed.trial.i)
  temp <- new.env()
  target.file <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%s.RData", seed)
  load(target.file, envir = temp)
  names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
  DDX1.pos.i <- which(grepl(pattern = "ENSG00000079785", x = names(temp[[sprintf("BMAseq.eFDR.Main.%s", data.type)]][["BMI"]])))
  DDX1.pos <- c(DDX1.pos, DDX1.pos.i)
  rm("temp")
}

names(DDX1.pos) <- seed.trial
DDX1.pos |> data.frame() |> View()

##------AZI2 position in the ranking list------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
best.model <- best.model |> column_to_rownames(var = "...1")
# SYMBOL         ENSEMBL
# 1   AZI2 ENSG00000163512

# General case for all possible seeds
# Find which random seed trial shows AZI2 as the cDEG detected with BMAseq at a 5000 ranking threshold
AZI2.best.model <- best.model[which(grepl(pattern = "ENSG00000163512", x = rownames(best.model))), 1:20]
AZI2.best.model <- AZI2.best.model |> unlist() # Transform the dataframe into a named vector
seed.trial <- character()
for (name in names(AZI2.best.model)) {
  if (!is.na(AZI2.best.model[name])) {
    seed.trial.i <- sub("bestmodel.", "", name)
    seed.trial <- c(seed.trial, seed.trial.i)
  }
}

var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
AZI2.pos <- numeric()
for (seed.trial.i in seed.trial) {
  seed <- sub(".*\\.", "", seed.trial.i)
  data.type <- sub("\\..*", "", seed.trial.i)
  temp <- new.env()
  target.file <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%s.RData", seed)
  load(target.file, envir = temp)
  names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
  AZI2.pos.i <- which(grepl(pattern = "ENSG00000163512", x = names(temp[[sprintf("BMAseq.eFDR.Main.%s", data.type)]][["BMI"]])))
  AZI2.pos <- c(AZI2.pos, AZI2.pos.i)
  rm("temp")
}

names(AZI2.pos) <- seed.trial
AZI2.pos |> data.frame() |> View()


##------CXCL10 position in the ranking list------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
best.model <- best.model |> column_to_rownames(var = "...1")
# SYMBOL         ENSEMBL
# 1   CXCL10 ENSG00000169245

# General case for all possible seeds
# Find which random seed trial shows CXCL10 as the cDEG detected with BMAseq at a 5000 ranking threshold
CXCL10.best.model <- best.model[which(grepl(pattern = "ENSG00000169245", x = rownames(best.model))), 1:20]
CXCL10.best.model <- CXCL10.best.model |> unlist() # Transform the dataframe into a named vector
seed.trial <- character()
for (name in names(CXCL10.best.model)) {
  if (!is.na(CXCL10.best.model[name])) {
    seed.trial.i <- sub("bestmodel.", "", name)
    seed.trial <- c(seed.trial, seed.trial.i)
  }
}

var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
CXCL10.pos <- numeric()
for (seed.trial.i in seed.trial) {
  seed <- sub(".*\\.", "", seed.trial.i)
  data.type <- sub("\\..*", "", seed.trial.i)
  temp <- new.env()
  target.file <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%s.RData", seed)
  load(target.file, envir = temp)
  names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
  CXCL10.pos.i <- which(grepl(pattern = "ENSG00000169245", x = names(temp[[sprintf("BMAseq.eFDR.Main.%s", data.type)]][["BMI"]])))
  CXCL10.pos <- c(CXCL10.pos, CXCL10.pos.i)
  rm("temp")
}

names(CXCL10.pos) <- seed.trial
CXCL10.pos |> data.frame() |> View()


uniq.cDEG <- c("IRF3", "DDX1", "AZI2", "CXCL10")
a <- lapply(uniq.cDEG, function(i) {name <- paste0(i, ".pos"); df <- get(name) |> data.frame(); df$ID <- rownames(df); colnames(df)[1] <- i; return(df)})
# CXCL10.pos <- CXCL10.pos |> data.frame() 
# colnames(CXCL10.pos) <- "CXCL10"

library(dplyr)
full_join(a[[1]], a[[2]], by = "ID") |> 
  full_join(a[[3]], by = "ID") |>
  full_join(a[[4]], by = "ID") |>
  dplyr::select(ID, everything()) |> View()
