## The 1st part - GOAL for one instance: build a BMAseq best model list per variable per random seed per ranking threshold
# Step 1: Load data from one seed
load("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti120FDR1.RData") 

# Step 2: Extract BMI-related index.DEG/name.DEG/best.model results using the training/test set in the random seed trial 120
temp.train <- output.multi.train.s2[["DEG.bestmodel"]][["BMI"]]
temp.test <- output.multi.test.s2[["DEG.bestmodel"]][["BMI"]]

# Step3: Identify BMI-related cDEGs among top 2000 ranked genes in both the training and test sets in the random seed trial 120
DEG.train <- BMAseq.eFDR.Main.train$BMI[1:2000] |> names() ## ENSG ID
DEG.test <- BMAseq.eFDR.Main.test$BMI[1:2000] |> names()
cDEG.all <- intersect(DEG.train, DEG.test) ## ENSG ID

# Step 4: Identify BMI-related unique DEGs among top 2000 ranked genes using the training/test set in the random seed trial 120
noncDEG.train <- setdiff(DEG.train, cDEG.all) ## ENSG ID
noncDEG.test <- setdiff(DEG.test, cDEG.all)

# Step 5: Extract index.DEG/name.DEG/best.model results mapped to each BMI-related cDEG among top 2000 ranked genes, using the training/test set in the random seed trial 120
cDEG.bestmodel.train <- temp.train[temp.train[, "name.DEG"] %in% cDEG.all, ] 
cDEG.bestmodel.test <- temp.test[temp.test[, "name.DEG"] %in% cDEG.all, ] 

# Step 6: Extract index.DEG/name.DEG/best.model results mapped to each BMI-related non-cDEG among top 2000 ranked genes, using the training/test set in the random seed trial 120
noncDEG.bestmodel.train <- temp.train[temp.train[, "name.DEG"] %in% noncDEG.train, ] 
noncDEG.bestmodel.test <- temp.test[temp.test[, "name.DEG"] %in% noncDEG.test, ] 

# Step 7: Put the best.model result in association with each BMI-related cDEG among top 2000 ranked genes using the training set in the random seed trial 120 and 
# the best.model result in association with each BMI-related cDEG among top 2000 ranked genes using the test set in the random seed trial 120 together
colnames(cDEG.bestmodel.train)[3] <- paste0("bestmodel.train", ".seed120") ## Change the name of the 3rd column that matches the element [cDEG.bestmodel.i] name of the large list [cDEG.bestmodel.list]; + line 166
colnames(cDEG.bestmodel.test)[3] <- paste0("bestmodel.test", ".seed120")
cDEG.bestmodel.all <- merge(cDEG.bestmodel.train, cDEG.bestmodel.test, by = c("index.DEG", "name.DEG")) ## merge: horizontal combination

# Step 8: Put the best.model result in association with each BMI-related non-cDEG among top 2000 ranked genes using the training set in the random seed trial 120 and 
# the best.model result in association with each BMI-related non-cDEG among top 2000 ranked genes using the test set in the random seed trial 120 together
colnames(noncDEG.bestmodel.train)[3] <- paste0("bestmodel", ".seed120")
colnames(noncDEG.bestmodel.test)[3] <- paste0("bestmodel", ".seed120")
noncDEG.bestmodel.all <- merge(noncDEG.bestmodel.train, noncDEG.bestmodel.test, all = T) 


## The 2nd part - GOAL for all instances: build a BMAseq best model list per variable per random seed per ranking threshold
# Step 1: Wrap up everything into a function with the following goals: 
# 1) build a BMAseq best model list for each cDEG in association with a variable among a ranking threshold using the training/test set in a random seed trial 
# 2) build a BMAseq best model list for each non-cDEG in association with a variable among a ranking threshold using the training/test set in a random seed trial 
# parameter: var.name [variable], threshold [ranking threshold], seed [random seed trial]
# input: parameter + data objects (output.multi.train.s2, output.multi.test.s2, BMAseq.eFDR.Main.train, BMAseq.eFDR.Main.test)
# output: 2 lists: 1) cDEG.model.result, 2) noncDEG.model.result 
collect_bestmodel <- function (var.name = NULL, threshold = NULL, seed = NULL) {
  temp.train <- output.multi.train.s2[["DEG.bestmodel"]][[var.name]]
  temp.test <- output.multi.test.s2[["DEG.bestmodel"]][[var.name]]
  DEG.train <- BMAseq.eFDR.Main.train[[var.name]][1:threshold] |> names()
  DEG.test <- BMAseq.eFDR.Main.test[[var.name]][1:threshold] |> names()
  cDEG.all <- intersect(DEG.train, DEG.test)
  noncDEG.train <- setdiff(DEG.train, cDEG.all)
  noncDEG.test <- setdiff(DEG.test, cDEG.all)
  cDEG.bestmodel.train <- temp.train[temp.train[, "name.DEG"] %in% cDEG.all, ] 
  cDEG.bestmodel.test <- temp.test[temp.test[, "name.DEG"] %in% cDEG.all, ] 
  noncDEG.bestmodel.train <- temp.train[temp.train[, "name.DEG"] %in% noncDEG.train, ] 
  noncDEG.bestmodel.test <- temp.test[temp.test[, "name.DEG"] %in% noncDEG.test, ] 
  colnames(cDEG.bestmodel.train)[3] <- sprintf("bestmodel.train.%s", seed)
  colnames(cDEG.bestmodel.test)[3] <- sprintf("bestmodel.test.%s", seed)
  colnames(noncDEG.bestmodel.train)[3] <- sprintf("bestmodel.train.%s", seed)
  colnames(noncDEG.bestmodel.test)[3] <- sprintf("bestmodel.test.%s", seed)
  cDEG.bestmodel.all <- merge(cDEG.bestmodel.train, cDEG.bestmodel.test, by = c("index.DEG", "name.DEG"))
  noncDEG.bestmodel.all <- merge(noncDEG.bestmodel.train, noncDEG.bestmodel.test, all = T)
  return(list(cDEG.model.result = cDEG.bestmodel.all, noncDEG.model.result = noncDEG.bestmodel.all))
}

# Step 2: Run the best model collector for each variable, each threshold, each seed
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
      cDEG.bestmodel.list <- c(cDEG.bestmodel.list, cDEG.bestmodel.i)
      noncDEG.bestmodel.list <- c(noncDEG.bestmodel.list, noncDEG.bestmodel.i)
    }
  }
}

# Step 3: save all best models of cDEGs and non-cDEGs among each variable, each threshold, each seed
date.analysis <- format(Sys.Date(), "%Y%b%d")
save(cDEG.bestmodel.list, noncDEG.bestmodel.list, file = sprintf("../ApplicationResult/Multi/RandomSeed/BestModelProportion/%s_BMAseq.bestmodel.list.RData", date.analysis))
load("../ApplicationResult/Multi/RandomSeed/BestModelProportion/2023May10_BMAseq.bestmodel.list.RData")


# The 3rd part - GOAL: 1) assign the most frequent BMAseq best model to each cDEG in association with a variable under a threshold, across 10 random seed trials
# 2) organize the results into a matrix:
# where each row represents a cDEG and the first 20 columns represent 10 pairs of training and test sets in 10 random seed trials
# the 21st column represents the most frequent BMAseq best model mapped to each cDEG
# the 22nd column represents the percentage of the most frequent BMAseq best model that is univariable among all models in the 21st column 
# the 23rd column represents the percentage of the most frequent BMAseq best model that is multivariable among all models in the 21st column
cDEG = cDEG.list = bestmodel.list = NULL
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- seq(1000, 5000, 1000)
most_fr_bestmodel_cDEG_per_var <- function(var.name = "BMI") {
  for (threshold.i in threshold.vec) {  
    for (seed.i in seed.vec) { 
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%sFDR1.RData", seed.i)
      load(file.name)
      # Free the memory
      object.list <- ls(pattern = "dat.*|.*s1")
      rm(list = object.list)
      # Extract variable-related index.DEG/name.DEG/best.model results using the training/test set in a random seed trial
      temp.train <- output.multi.train.s2[["DEG.bestmodel"]][[var.name]]
      temp.test <- output.multi.test.s2[["DEG.bestmodel"]][[var.name]]
      cDEG <- intersect(BMAseq.eFDR.Main.train[[var.name]][1:threshold.i] |> names(), BMAseq.eFDR.Main.test[[var.name]][1:threshold.i] |> names())
      cDEG.list <- c(cDEG.list, list(cDEG))
    }
    cDEG.unique <- unique(unlist(cDEG.list)) 
    
    # Create a matrix that shows the best model of each unique gene among the 10 random seed trials, in a threshold 
    cDEG.bestmodel.matrix <- matrix(0, nrow = length(cDEG.unique), ncol = 20) 
    ## nrow: each row represents a unique cDEG across 10 random seed trials, in a threshold
    ## ncol: 10 pairs of training and test sets in 10 random seed trials, in a threshold 
    rownames(cDEG.bestmodel.matrix) <- cDEG.unique
    colnames(cDEG.bestmodel.matrix) <- paste(c("bestmodel.train", "bestmodel.test"), rep(seed.vec, each = 2), sep = ".")  
    ## "bestmodel.train.8809678" "bestmodel.test.8809678" "bestmodel.train.98907" "bestmodel.test.98907" ...
    for (i in 1:10) { # Loop through 10 seeds
      seed.value <- seed.vec[i]
      temp.cDEG.bestmodel <- cDEG.bestmodel.list[[sprintf("%s+%s+%s", var.name, threshold.i, seed.value)]]
      ## `cDEG.bestmodel.list` contains the best models of cDEGs among all variables, seeds, thresholds 
      ## Each element in this list contains the best models of cDEGs in association with a variable in a threshold, using the training and test sets in a random seed trial
      ## For the purpose of extraction, the ID for each element is composed of var.name, threshold.i, and seed.value
      ## temp.cDEG.bestmodel == best models of cDEGs in association with a variable in a threshold, using the training and test sets in a random seed trial
      
      temp.noncDEG.bestmodel <- noncDEG.bestmodel.list[[sprintf("%s+%s+%s", var.name, threshold.i, seed.value)]]
      ## `noncDEG.bestmodel.list` contains the best models of non-cDEGs among all variables, seeds, thresholds 
      ## Each element in this list contains the best models of non-cDEGs in association with a variable in a threshold, using the training and test sets in a random seed trial
      ## For the purpose of extraction, the ID for each element is composed of var.name, threshold.i, and seed.value
      ## temp.noncDEG.bestmodel == best models of non-cDEGs in association with a variable in a threshold, using the training and test sets in a random seed trial
      
      # Extract the ENSG IDs of cDEGs/non-cDEGs in association with a variable in a threshold, using the training and test sets in a random seed trial
      cDEG.per.seed <- temp.cDEG.bestmodel$name.DEG ## ENSG ID
      noncDEG.per.seed <- temp.noncDEG.bestmodel$name.DEG ## ENSG ID
      
      for (j in 1:length(cDEG.unique)) { ## cDEG.unique: unique cDEGs across 10 random seed trials, in a threshold
        if (cDEG.unique[j] %in% cDEG.per.seed) { ## cDEG.unique[j] is the top ranked cDEG in this seed trial, under a particular threshold
          cDEG.bestmodel.matrix[j, 2*i-1] <- temp.cDEG.bestmodel[which(temp.cDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.train.%s", seed.value)] # Column name maps to a seed value
          ## cDEG.bestmodel.matrix <- matrix(0, nrow = length(cDEG.unique), ncol = 20)
          ## 2*i - 1: training; 2*i: test
          ## each row represents a unique cDEG across 10 random seed trials, in a threshold 
          ## temp.cDEG.bestmodel == best models of cDEGs in association with a variable in a threshold, using the training and test sets in a random seed trial
          cDEG.bestmodel.matrix[j, 2*i] <- temp.cDEG.bestmodel[which(temp.cDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.test.%s", seed.val)]
        } else if (cDEG.unique[j] %in% noncDEG.per.seed) { ## cDEG.unique[j] is the top ranked non-cDEG in this seed trial
          cDEG.bestmodel.matrix[j, 2*i-1] <- temp.noncDEG.bestmodel[which(temp.noncDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.train.%s", seed.val)]
          cDEG.bestmodel.matrix[j, 2*i] <- temp.noncDEG.bestmodel[which(temp.noncDEG.bestmodel$name.DEG == cDEG.unique[j]), sprintf("bestmodel.test.%s", seed.val)]
        } else { ## cDEG.unique[j] is neither the top ranked cDEG nor the top ranked non-cDEG in this seed trial, under a particular threshold
          ## in other words, cDEG.unique[j] is not a top ranked DEG in this seed trial, under a particular threshold
          ## 2023Apr29_BMAseq_ManuscriptResult_ALiu_iMac slice 18: number of age-related cDEGs under the ranking threshold 4000 does not include the age-related cDEGs under the other ranking threshold
          cDEG.bestmodel.matrix[j, 2*i-1] <- NA
          cDEG.bestmodel.matrix[j, 2*i] <- NA
        }
      }
    }
    
    # Return a list where each element is the most frequent best model formula for each cDEG
    most.fr.bestmodel <- apply(cDEG.bestmodel.matrix, 1, function(x) names(table(x))[which(table(x) == max(table(x)))]) 
    
    # If there are multiple best models, e.g., univariable + multivariable; multivariable + multivariable
    most.fr.bestmodel.onecol <- lapply(most.fr.bestmodel, function(x) paste0(x, collapse = "/")) |> unlist()
    
    # Sum of most frequent univariable best models across all cDEGs
    sum.uni <- lapply(most.fr.bestmodel, function (x) any(sprintf("~1+%s", var.name) %in% x)) |> unlist() |> sum()
    
    # If the best model is multivariable [x != sprintf("~1+%s", var.name] and contains the target phenotype [grepl(var.name, x)]
    bestmodel.multi <- lapply(most.fr.bestmodel, function (x) { (x != sprintf("~1+%s", var.name)) & grepl(var.name, x) }) # lapply returns a list; each element contains the logical object (TRUE/FALSE)
    
    # If there are multiple multivariable best models for a cDEG 
    sum.multi <- lapply(bestmodel.multi, sum)
    ## lapply returns a list; each element is the sum of most frequent multivariable best models per cDEG
    
    # Sum of most frequent multivariable best models across all cDEGs
    sum.multi <- ifelse(sum.multi > 1, 1, sum.multi) |> unlist() |> sum() 
    
    # Percentage of most frequent univariable best models across all the most frequent best models containing the target variable
    percent.uni <- sum.uni / sum(sum.uni + sum.multi)
    
    # Percentage of most frequent multivariable best models across all the most frequent best models containing the target variable
    percent.multi <- sum.multi / sum(sum.uni + sum.multi)
    
    cDEG.bestmodel.matrix <- cbind(cDEG.bestmodel.matrix, most.fr.bestmodel.onecol, sum.uni, sum.multi, percent.uni, percent.multi)
    file.name <- sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.RDS", var.name, threshold.i)
    saveRDS(cDEG.bestmodel.matrix, file.name)
    print(paste("Complete the BMAseq best model matrix of cDEG associated with", var.name, "among 10 random seeds for the threshold", threshold.i))
    cDEG = cDEG.list =  NULL
    
    # Free the memory
    object.list <- ls(pattern = "BMAseq.*|.*s2")
    rm(list = object.list)
  }
}  

mapply(FUN = most_fr_bestmodel_cDEG_per_var,
       c("BMI", "AGE", "SEX", "MHABNWBC"))

# Spot check the data 
BMI_1000 <- readRDS("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_1000.RDS")
temp.bestmodel <- BMAseq.bestmodel.list[["BMI+1000+8809678"]]


# The 4th part - GOAL: export each matrix which captures the most frequent BMAseq best model to each cDEG in association with a variable under a threshold, across 10 random seed trials
# Export the matrix into xslx
library(openxlsx)
threshold.vec <- seq(1000, 5000, 1000)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
for (threshold.i in threshold.vec) {
  for (var.i in var.vec) {
    file.name <- sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.RDS", var.i, threshold.i)
    data <- readRDS(file.name) |> data.frame()
    file.name <- sprintf("../ApplicationData/derived/RandomSeed/BestModelMatrix/%s_%s.xlsx", var.i, threshold.i)
    write.xlsx(x = data, file = file.name, rowNames = T, colWidths = 40, firstActiveRow = 2, firstActiveCol = 2)
  }
}
