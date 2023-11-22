# Load data files
library(readxl)
data_dir <- '../data/'
file_names <- c(
  paste0(data_dir, 'pheno_data.xlsx'),
  paste0(data_dir, 'dict_data1.xlsx'),
  paste0(data_dir, 'dict_data2.xlsx')
)
loaded_data <- lapply(file_names, read_xlsx) # Load xlsx files into a list
names(loaded_data) <- c('pheno_data', 'dict_data1', 'dict_data2')
saveRDS(loaded_data, file=paste0(data_dir, 'loaded_data.RDS'))

loaded_data <- readRDS(paste0(data_dir, 'loaded_data.RDS'))
pheno_data <- loaded_data$pheno_data
dict_data1 <- loaded_data$dict_data1
dict_data2 <- loaded_data$dict_data2
dict_full <- rbind(dict_data1, dict_data2)
names(pheno_data)
table(pheno_data$SMTSD, useNA='ifany') # View the distribution of tissue type
idx_sub_fat <- which(pheno_data$SMTSD=='Adipose - Subcutaneous')
dict_name <- dict_full$name
length(dict_name); ncol(pheno_data)
