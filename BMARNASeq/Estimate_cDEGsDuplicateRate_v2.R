cDEGs = cDEGs.list = NULL
threshold.vec <- seq(1000, 5000, 1000)
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)

duplicate_rate_cDEGs_per_variable <- function(var.name = "BMI", method.name = "DESeq2", model.type = "Uni") {
  for (threshold.i in threshold.vec) {
    for (seed.i in seed.vec) {
      if (method.name == "BMAseq" & var.name != "BMIxSEX") {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%s%s%s.RData", method.name, model.type, seed.i))
        cDEGs <- c(intersect(names(get(sprintf("%s.eFDR.Main.train", method.name))[[var.name]][1:threshold.i]), 
                             names(get(sprintf("%s.eFDR.Main.test", method.name))[[var.name]][1:threshold.i])))
      } else if (method.name == "BMAseq" & var.name == "BMIxSEX") {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%s%s%s.RData", method.name, model.type, seed.i))
        var.name2 <- "BMI"
        cDEGs <- c(intersect(names(get(sprintf("%s.eFDR.Interaction.train", method.name))[[var.name2]][1:threshold.i]), 
                             names(get(sprintf("%s.eFDR.Interaction.test", method.name))[[var.name2]][1:threshold.i])))
      } else if (method.name %in% c("DESeq2", "edgeR")) {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%s%s%s.RData", method.name, model.type, seed.i))
        cDEGs <- c(intersect(get(sprintf("%s.eFDR.GeneName.train", method.name))[[var.name]][1:threshold.i], 
                             get(sprintf("%s.eFDR.GeneName.test", method.name))[[var.name]][1:threshold.i]))
      } else if (method.name == "eBayes") {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%s%s%s.RData", method.name, model.type, seed.i))
        cDEGs <- c(intersect(names(get(sprintf("%s.eFDR.train2", method.name))[[var.name]][1:threshold.i]), 
                             names(get(sprintf("%s.eFDR.test2", method.name))[[var.name]][1:threshold.i])))
      } else {
        load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/%s%s%s.RData", method.name, model.type, seed.i))
        method.name2 <- "voom"
        cDEGs <- c(intersect(names(get(sprintf("%s.eFDR.train2", method.name2))[[var.name]][1:threshold.i]), 
                             names(get(sprintf("%s.eFDR.test2", method.name2))[[var.name]][1:threshold.i])))
      }
      cDEGs.list <- append(cDEGs.list, list(cDEGs))
    }
    cDEGs.unique <- unique(unlist(cDEGs.list))
    
    # Create a matrix that shows the occurrence of each unique gene among the 10 random seed trials
    cDEGs.matrix <- matrix(0, nrow = length(cDEGs.unique), ncol = 10)
    rownames(cDEGs.matrix) <- cDEGs.unique
    colnames(cDEGs.matrix) <- seed.vec   
    for (i in 1:10) {
      cDEGs.i.seed <- cDEGs.list[[i]]
      for (j in 1:length(cDEGs.unique)) {
        if (cDEGs.unique[j] %in% cDEGs.i.seed) {
          cDEGs.matrix[j, i] <- 1
        }
      }
    }
    Sum <- apply(cDEGs.matrix, 1, sum)
    Percent <- Sum/10
    cDEGs.matrix <- cbind(cDEGs.matrix, Sum, Percent) # Add two new columns to the matrix
    saveRDS(cDEGs.matrix, sprintf("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/%s%s%s_%s.RDS", method.name, model.type, var.name, threshold.i))
    
    print(paste("Complete the duplicated rate matrix of cDEGs associated with", var.name, "identified by", method.name, "among 10 random seeds for the threshold", threshold.i))
    cDEGs = cDEGs.list =  NULL
  }
}

mapply(duplicate_rate_cDEGs_per_variable,
       rep(c("BMI", "AGE", "SEX", "MHABNWBC"), each = 4), 
       rep(c("DESeq2", "edgeR", "eBayes", "VoomLimma"), times = 4))


# For illustration purpose, only consider the phenotype BMI and 2000 ranking threshold
method.vec <- c("DESeq2", "edgeR", "eBayes", "VoomLimma")
model.type <- "Uni"
var.name <- "BMI"
threshold.i <- 2000
for (method.name in method.vec) {
  data <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/%s%s%s_%s.RDS", method.name, model.type, var.name, threshold.i))
  write.csv(data, sprintf("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/%s%s%s_%s.csv", method.name, model.type, var.name, threshold.i))
}
