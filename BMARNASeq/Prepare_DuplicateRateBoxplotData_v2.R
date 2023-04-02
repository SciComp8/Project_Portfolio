# Update on Prepare_DuplicateRateBoxplotData_v1.R: 
# To make the function `boxplot_duplicate_rate_data()` compatible with all modeling scenariors (univariable, multivariable, and interaction across all approaches),
# an argument named `model.type` is added.

boxplot_duplicate_rate_data <- function (var.name = NULL, method.name = NULL, model.type = NULL, threshold = NULL) {
  load.data <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/%s%s%s_%s.RDS", method.name, model.type, var.name, threshold))
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
