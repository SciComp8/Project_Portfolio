## The first part
# Load data from one seed
load("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti120FDR1.RData") 

# Extract index.DEG/name.DEG/best.model results related to BMI from the training set and test set
temp.train <- output.multi.train.s2[["DEG.bestmodel"]][["BMI"]]
temp.test <- output.multi.test.s2[["DEG.bestmodel"]][["BMI"]]

# Extract the BMI-related cDEGs among top 2000 ranked genes from the training set and among top 2000 ranked genes from the test set
cDEG.all <- intersect(BMAseq.eFDR.Main.train$BMI[1:2000] |> names(), BMAseq.eFDR.Main.test$BMI[1:2000] |> names())  # ranking threshold: 2000

# Extract index.DEG/name.DEG/best.model results mapped to BMI-related cDEGs from the training set and test set
cDEG.bestmodel.train <- temp.train[temp.train[, 2] %in% cDEG.all, ] # 2: name.DEG; save this object
cDEG.bestmodel.test <- temp.test[temp.test[, 2] %in% cDEG.all, ] 

# Combine cDEG best model results from the training set and test set
colnames(cDEG.bestmodel.train)[3] <- paste0("bestmodel.train", ".seed120")
colnames(cDEG.bestmodel.test)[3] <- paste0("bestmodel.test", ".seed120")
cDEG.bestmodel.all <- merge(cDEG.bestmodel.train, cDEG.bestmodel.test, by = c("index.DEG", "name.DEG"))

# Put everything into a function
collect_bestmodel <- function (threshold = NULL, seed = NULL, var.name = NULL) {
  temp.train <- output.multi.train.s2[["DEG.bestmodel"]][[var.name]]
  temp.test <- output.multi.test.s2[["DEG.bestmodel"]][[var.name]]
  cDEG.all <- intersect(BMAseq.eFDR.Main.train[[var.name]][1:threshold] |> names(), BMAseq.eFDR.Main.test[[var.name]][1:threshold] |> names())
  cDEG.bestmodel.train <- temp.train[temp.train[, 2] %in% cDEG.all, ] 
  cDEG.bestmodel.test <- temp.test[temp.test[, 2] %in% cDEG.all, ] 
  colnames(cDEG.bestmodel.train)[3] <- sprintf("bestmodel.train.%s", seed)
  colnames(cDEG.bestmodel.test)[3] <- sprintf("bestmodel.test.%s", seed)
  cDEG.bestmodel.all <- merge(cDEG.bestmodel.train, cDEG.bestmodel.test, by = c("index.DEG", "name.DEG"))
  return(list(model.result = cDEG.bestmodel.all))
}

bestmodel.list <- NULL
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- seq(1000, 5000, 1000)
for (threshold.i in threshold.vec) {
  for (seed.i in seed.vec) {
    load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%sFDR1.RData", seed.i))
    for (var.i in var.vec) {
      bestmodel.i <- collect_bestmodel(threshold = threshold.i, seed = seed.i, var.name = var.i)[["model.result"]] |> list()
      names(bestmodel.i) <- paste0(var.i, "+", threshold.i, "+", seed.i)
      bestmodel.list <- append(bestmodel.list, bestmodel.i)
    }
  }
}

BMI.bestmodel.1000 <- bestmodel.list[names(bestmodel.list) %in% "BMI+1000"]
date.analysis <- format(Sys.Date(), "%Y%b%d")
saveRDS(bestmodel.list, file = sprintf("../ApplicationResult/Multi/RandomSeed/BestModelProportion/%s_BMAseq.bestmodel.list.RDS", date.analysis))
BMAseq.bestmodel.list <- readRDS("../ApplicationResult/Multi/RandomSeed/BestModelProportion/2023May07_BMAseq.bestmodel.list.RDS")



## The second part
cDEGs = cDEGs.list = bestmodel.list = NULL
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- seq(1000, 5000, 1000)
bestmodel_cDEGs_per_variable <- function(var.name = "BMI") {
  for (threshold.i in threshold.vec) {
    for (seed.i in seed.vec) {
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%sFDR1.RData", seed.i))
      temp.train <- output.multi.train.s2[["DEG.bestmodel"]][[var.name]]
      temp.test <- output.multi.test.s2[["DEG.bestmodel"]][[var.name]]
      cDEGs <- intersect(BMAseq.eFDR.Main.train[[var.name]][1:threshold.i] |> names(), BMAseq.eFDR.Main.test[[var.name]][1:threshold.i] |> names())
      cDEGs.list <- append(cDEGs.list, list(cDEGs))
    }
    cDEGs.unique <- unique(unlist(cDEGs.list))
    
    # Create a matrix that shows the best model of each unique gene among the 10 random seed trials
    cDEGs.bestmodel.matrix <- matrix(0, nrow = length(cDEGs.unique), ncol = 20)
    rownames(cDEGs.bestmodel.matrix) <- cDEGs.unique
    colnames(cDEGs.bestmodel.matrix) <- paste(c("bestmodel.train", "bestmodel.test"), rep(seed.vec, each = 2))  
    for (i in 1:10) {
      cDEGs.i.seed <- cDEGs.list[[i]]
      i.seed <- seed.vec[i]
      temp.bestmodel <- BMAseq.bestmodel.list[[sprintf("%s+%s+%s", var.name, threshold.i, i.seed)]]
      for (j in 1:length(cDEGs.unique)) {
        if (cDEGs.unique[j] %in% cDEGs.i.seed) {
          cDEGs.bestmodel.matrix[j, 2*i-1] <- temp.bestmodel[which(temp.bestmodel$name.DEG == cDEGs.unique[j]), sprintf("bestmodel.train.%s", i.seed)]
          cDEGs.bestmodel.matrix[j, 2*i] <- temp.bestmodel[which(temp.bestmodel$name.DEG == cDEGs.unique[j]), sprintf("bestmodel.test.%s", i.seed)]
        }
      }
    }
    Sum.uni <- apply(cDEGs.bestmodel.matrix, 1, function(x) sum(x == sprintf("~1+%s", var.name)))
    Percent.uni <- Sum.uni/20
    Sum.multi <- apply(cDEGs.bestmodel.matrix, 1, function(x) {
      cDEG.bestmodel.split <- strsplit(x, split = "+", fixed = T)
      cDEG.bestmodel.split <- lapply(cDEG.bestmodel.split, function(x) x[-1])
      cDEG.bestmodel.split.filter <- cDEG.bestmodel.split[grepl(var.name, cDEG.bestmodel.split)] # Multivariable model formula contains the phenotype
      var.num <- unlist(lapply(cDEG.bestmodel.split.filter, function(x) length(x)))
      sum(var.num > 1)
    })
    Percent.multi <- Sum.multi/20
    cDEGs.bestmodel.matrix <- cbind(cDEGs.bestmodel.matrix, Sum.uni, Percent.uni, Sum.multi, Percent.multi)
    saveRDS(cDEGs.bestmodel.matrix, sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.RDS", var.name, threshold.i))
    print(paste("Complete the BMAseq best model matrix of cDEGs associated with", var.name, "among 10 random seeds for the threshold", threshold.i))
    cDEGs = cDEGs.list =  NULL
  }
}  

mapply(FUN = bestmodel_cDEGs_per_variable,
       c("BMI", "AGE", "SEX", "MHABNWBC"))

# Spot check the data 
BMI_1000 <- readRDS("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_1000.RDS")
temp.bestmodel <- BMAseq.bestmodel.list[["BMI+1000+98907"]]
