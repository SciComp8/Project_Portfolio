## The first part - build a BMAseq best model list per variable per random seed per ranking threshold
# Load data from one seed
load("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti120FDR1.RData") 

# Extract index.DEG/name.DEG/best.model results related to BMI from the training set and test set
temp.train <- output.multi.train.s2[["DEG.bestmodel"]][["BMI"]]
temp.test <- output.multi.test.s2[["DEG.bestmodel"]][["BMI"]]

# Revise: Extract the BMI-related cDEGs among top 2000 ranked genes from the training set and among top 2000 ranked genes from the test set
cDEG.train <- BMAseq.eFDR.Main.train$BMI[1:2000] |> names()
cDEG.test <- BMAseq.eFDR.Main.test$BMI[1:2000] |> names()
cDEG.all <- intersect(cDEG.train, cDEG.test)
noncDEG.train <- setdiff(cDEG.train, cDEG.all)
noncDEG.test <- setdiff(cDEG.test, cDEG.all)

# Extract index.DEG/name.DEG/best.model results mapped to BMI-related cDEGs from the training set and test set
cDEG.bestmodel.train <- temp.train[temp.train[, 2] %in% cDEG.all, ] # 2: name.DEG; save this object
cDEG.bestmodel.test <- temp.test[temp.test[, 2] %in% cDEG.all, ] 

# Combine cDEG best model results from the training set and test set
colnames(cDEG.bestmodel.train)[3] <- paste0("bestmodel.train", ".seed120")
colnames(cDEG.bestmodel.test)[3] <- paste0("bestmodel.test", ".seed120")
cDEG.bestmodel.all <- merge(cDEG.bestmodel.train, cDEG.bestmodel.test, by = c("index.DEG", "name.DEG"))

# Add: Extract index.DEG/name.DEG/best.model results mapped to BMI-related non-cDEGs from the training set and test set
noncDEG.bestmodel.train <- temp.train[temp.train[, 2] %in% noncDEG.train, ] # 2: name.DEG; save this object
noncDEG.bestmodel.test <- temp.test[temp.test[, 2] %in% noncDEG.test, ] 

# Add: Combine non-cDEG best model results from the training set and test set
colnames(noncDEG.bestmodel.train)[3] <- paste0("bestmodel", ".seed120")
colnames(noncDEG.bestmodel.test)[3] <- paste0("bestmodel", ".seed120")
noncDEG.bestmodel.all <- rbind(noncDEG.bestmodel.train, noncDEG.bestmodel.test)

# Put everything into a function
collect_bestmodel <- function (threshold = NULL, seed = NULL, var.name = NULL) {
  temp.train <- output.multi.train.s2[["DEG.bestmodel"]][[var.name]]
  temp.test <- output.multi.test.s2[["DEG.bestmodel"]][[var.name]]
  cDEG.train <- BMAseq.eFDR.Main.train[[var.name]][1:threshold] |> names()
  cDEG.test <- BMAseq.eFDR.Main.test[[var.name]][1:threshold] |> names()
  cDEG.all <- intersect(cDEG.train, cDEG.test)
  noncDEG.train <- setdiff(cDEG.train, cDEG.all)
  noncDEG.test <- setdiff(cDEG.test, cDEG.all)
  cDEG.bestmodel.train <- temp.train[temp.train[, 2] %in% cDEG.all, ] 
  cDEG.bestmodel.test <- temp.test[temp.test[, 2] %in% cDEG.all, ] 
  colnames(cDEG.bestmodel.train)[3] <- sprintf("bestmodel.train.%s", seed)
  colnames(cDEG.bestmodel.test)[3] <- sprintf("bestmodel.test.%s", seed)
  cDEG.bestmodel.all <- merge(cDEG.bestmodel.train, cDEG.bestmodel.test, by = c("index.DEG", "name.DEG"))
  noncDEG.bestmodel.train <- temp.train[temp.train[, 2] %in% noncDEG.train, ] # 2: name.DEG; save this object
  noncDEG.bestmodel.test <- temp.test[temp.test[, 2] %in% noncDEG.test, ] 
  colnames(noncDEG.bestmodel.train)[3] <- sprintf("bestmodel.train.%s", seed)
  colnames(noncDEG.bestmodel.test)[3] <- sprintf("bestmodel.test.%s", seed)
  noncDEG.bestmodel.all <- merge(noncDEG.bestmodel.train, noncDEG.bestmodel.test, all = T)
  return(list(cDEG.model.result = cDEG.bestmodel.all, noncDEG.model.result = noncDEG.bestmodel.all))
}

cDEG.bestmodel.list = noncDEG.bestmodel.list = NULL
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- seq(1000, 5000, 1000)
for (threshold.i in threshold.vec) {
  for (seed.i in seed.vec) {
    load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%sFDR1.RData", seed.i))
    for (var.i in var.vec) {
      cDEG.bestmodel.i <- collect_bestmodel(threshold = threshold.i, seed = seed.i, var.name = var.i)[["cDEG.model.result"]] |> list()
      noncDEG.bestmodel.i <- collect_bestmodel(threshold = threshold.i, seed = seed.i, var.name = var.i)[["noncDEG.model.result"]] |> list()
      names(cDEG.bestmodel.i) = names(noncDEG.bestmodel.i) = paste0(var.i, "+", threshold.i, "+", seed.i)
      cDEG.bestmodel.list <- append(cDEG.bestmodel.list, cDEG.bestmodel.i)
      noncDEG.bestmodel.list <- append(noncDEG.bestmodel.list, noncDEG.bestmodel.i)
    }
  }
}

date.analysis <- format(Sys.Date(), "%Y%b%d")
save(cDEG.bestmodel.list, noncDEG.bestmodel.list, file = sprintf("../ApplicationResult/Multi/RandomSeed/BestModelProportion/%s_BMAseq.bestmodel.list.RData", date.analysis))
load("../ApplicationResult/Multi/RandomSeed/BestModelProportion/2023May10_BMAseq.bestmodel.list.RData")

## The second part - assign each BMAseq best model to each cDEG per variable per random seed per ranking threshold
cDEG = cDEG.list = bestmodel.list = NULL
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- seq(1000, 5000, 1000)
most_fr_bestmodel_cDEG_per_var <- function(var.name = "BMI") {
  for (threshold.i in threshold.vec) {
    for (seed.i in seed.vec) {
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%sFDR1.RData", seed.i))
      temp.train <- output.multi.train.s2[["DEG.bestmodel"]][[var.name]]
      temp.test <- output.multi.test.s2[["DEG.bestmodel"]][[var.name]]
      cDEG <- intersect(BMAseq.eFDR.Main.train[[var.name]][1:threshold.i] |> names(), BMAseq.eFDR.Main.test[[var.name]][1:threshold.i] |> names())
      cDEG.list <- append(cDEG.list, list(cDEG))
    }
    cDEG.unique <- unique(unlist(cDEG.list))
    
    # Create a matrix that shows the best model of each unique gene among the 10 random seed trials
    cDEG.bestmodel.matrix <- matrix(0, nrow = length(cDEG.unique), ncol = 20)
    rownames(cDEG.bestmodel.matrix) <- cDEG.unique
    colnames(cDEG.bestmodel.matrix) <- paste(c("bestmodel.train", "bestmodel.test"), rep(seed.vec, each = 2))  
    for (i in 1:10) { # Loop through seed
      seed.val <- seed.vec[i]
      temp.cDEG.bestmodel <- cDEG.bestmodel.list[[sprintf("%s+%s+%s", var.name, threshold.i, seed.val)]]
      temp.noncDEG.bestmodel <- noncDEG.bestmodel.list[[sprintf("%s+%s+%s", var.name, threshold.i, seed.val)]]
      cDEG.per.seed <- temp.cDEG.bestmodel$name.DEG
      noncDEG.per.seed <- temp.noncDEG.bestmodel$name.DEG
      for (j in 1:length(cDEG.unique)) {
        if (cDEG.unique[j] %in% cDEG.per.seed) {
          cDEG.bestmodel.matrix[j, 2*i-1] <- temp.cDEG.bestmodel[which(temp.cDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.train.%s", seed.val)] # Column name
          cDEG.bestmodel.matrix[j, 2*i] <- temp.cDEG.bestmodel[which(temp.cDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.test.%s", seed.val)]
        } else if (cDEG.unique[j] %in% noncDEG.per.seed) {
          cDEG.bestmodel.matrix[j, 2*i-1] <- temp.noncDEG.bestmodel[which(temp.noncDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.train.%s", seed.val)]
          cDEG.bestmodel.matrix[j, 2*i] <- temp.noncDEG.bestmodel[which(temp.noncDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.test.%s", seed.val)]
        } else { # DEG is not the top ranked DEG in this seed trial
          cDEG.bestmodel.matrix[j, 2*i-1] <- NA
          cDEG.bestmodel.matrix[j, 2*i] <- NA
        }
      }
    }
    most.fr.bestmodel <- apply(cDEG.bestmodel.matrix, 1, function(x) names(table(x))[which(table(x) == max(table(x)))]) # Return a list where each element is the most frequent best model formula
    most.fr.bestmodel.onecol <- lapply(most.fr.bestmodel, function(x) paste0(x, collapse = "/")) |> unlist() # If there are multiple best models
    sum.uni <- lapply(most.fr.bestmodel, function (x) any(sprintf("~1+%s", var.name) %in% x)) |> unlist() |> sum()
    percent.uni <- sum.uni/length(most.fr.bestmodel)
    bestmodel.multi <- lapply(most.fr.bestmodel, function (x) { grepl(var.name, x) & (x != sprintf("~1+%s", var.name)) }) 
    sum.multi <- lapply(bestmodel.multi, sum) # If there are multiple multivariable best models ... 
    sum.multi <- ifelse(sum.multi > 1, 1, sum.multi) |> unlist() |> sum()
    percent.multi <- sum.multi/length(most.fr.bestmodel)
    cDEG.bestmodel.matrix <- cbind(cDEG.bestmodel.matrix, most.fr.bestmodel.onecol, percent.uni, percent.multi)
    saveRDS(cDEG.bestmodel.matrix, sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.RDS", var.name, threshold.i))
    print(paste("Complete the BMAseq best model matrix of cDEG associated with", var.name, "among 10 random seeds for the threshold", threshold.i))
    cDEG = cDEG.list =  NULL
  }
}  

mapply(FUN = most_fr_bestmodel_cDEG_per_var,
       c("BMI", "AGE", "SEX", "MHABNWBC"))

# Spot check the data 
BMI_1000 <- readRDS("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_1000.RDS")
temp.bestmodel <- BMAseq.bestmodel.list[["BMI+1000+8809678"]]

# Export the results into xslx
library(openxlsx)
threshold.vec <- seq(1000, 5000, 1000)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
for (threshold.i in threshold.vec) {
  for (var.i in var.vec) {
    data <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.RDS", var.i, threshold.i)) |> data.frame()
    write.xlsx(data, sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.xlsx", var.i, threshold.i), rowNames = T, colWidths = 40, firstActiveRow = 2, firstActiveCol = 2)
  }
}
