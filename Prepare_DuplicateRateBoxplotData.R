library(data.table)
boxplot_duplicate_rate_data <- function (var.name = NULL, method.name = NULL, threshold = NULL) {
  load.data <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/%sMultiInt%s_%s.RDS", method.name, var.name, threshold))
  temp.data <- data.table(threshold=threshold, 
                          method=method.name,
                          variable=var.name,
                          cDEGs.duplicate.rate.median=median(load.data[,12]),
                          cDEGs.duplicate.rate.mean=mean(load.data[,12]))
  return(temp.data)
}

organized_duplicate_rate_data <- mapply(boxplot_duplicate_rate_data,
                                        var.name = rep(c("BMI", "AGE","SEX", "MHABNWBC", "BMIxSEX"), each = 25), 
                                        method.name = rep(rep(c("BMAseq", "DESeq2", "edgeR", "eBayes", "VoomLimma"), each = 5), times = 5),
                                        threshold = rep(rep(seq(1000, 5000, 1000), times = 5), times = 5)) |> t()
rownames(organized_duplicate_rate_data) <- NULL
organized_duplicate_rate_data <- apply(organized_duplicate_rate_data, 2, unlist)
organized_duplicate_rate_data <- as.data.table(organized_duplicate_rate_data)
