# Load data files
library(readxl)
library(tidyverse)

# >
data_dir <- '../data/'
file_names <- c(
  paste0(data_dir, 'pheno_data.xlsx'),
  paste0(data_dir, 'dict_data1.xlsx'),
  paste0(data_dir, 'dict_data2.xlsx')
)
loaded_data <- lapply(file_names, read_xlsx, na = c('NA')) # Load xlsx files into a list
names(loaded_data) <- c('pheno_data', 'dict_data1', 'dict_data2')
saveRDS(loaded_data, file=paste0(data_dir, 'loaded_data.RDS'))

# >
loaded_data <- readRDS(paste0(data_dir, 'loaded_data.RDS'))
pheno_data <- loaded_data$pheno_data
dict_data1 <- loaded_data$dict_data1
dict_data2 <- loaded_data$dict_data2
dict_full <- rbind(dict_data1, dict_data2)

# Clean data
names(pheno_data)
table(pheno_data$SMTSD, useNA='ifany') # View the distribution of tissue type
idx_sub_fat <- which(pheno_data$SMTSD=='Adipose - Subcutaneous')
dict_name <- dict_full$name
length(dict_name); ncol(pheno_data)
colnames(pheno_data)[!colnames(pheno_data) %in% dict_name] # [1] "dbGaP_Subject_ID.x" "dbGaP_Subject_ID.y" "dbGaP_Sample_ID" 
pheno_data_filter <- pheno_data[, colnames(pheno_data) %in% dict_name]
var_str <- dict_full$type=='string'
sum(var_str) # Number of string variables
var_focus <- c('COHORT', as.character(dict_full$name[!var_str]))
var_cat <- c(1, as.numeric(dict_full$type[!var_str] == 'integer, encoded value'))
colSums(is.na(pheno_data_filter[, var_focus]))
var_remove <- which(colSums(is.na(pheno_data_filter[, var_focus])) > 70)
var_focus_2 <- var_focus[-var_remove]; var_cat_2 <- var_cat[-var_remove]
# duplicated(c(2, 2, 2))
# [1] FALSE  TRUE  TRUE

# Identify columns in the pheno_data_filter dataframe where all elements are repetitive 
row_count_minus_one <- nrow(pheno_data_filter) - 1
idx_identity <- vector()
for (col in var_focus_2) {
  sum_dup <- sum(duplicated(pheno_data_filter[, col])) # Calculate the sum of duplicated elements in the column
  if (sum_dup == row_count_minus_one) { # Check if the sum of duplicates is equal to the threshold
    idx_identity <- c(idx_identity, which(names(pheno_data_filter) == col)) # Store the index if condition is met
  }
}
idx_identity # [1] 222 225 227 228 229 230 233 239 242 243 244 261 262
var_focus_3 <- var_focus_2[-idx_identity]; var_cat_3 <- var_cat_2[-idx_identity]
taget_col <- colnames(pheno_data_filter[var_focus_3[!as.logical(var_cat_3)]])
pheno_data_filter <- pheno_data_filter |> 
  mutate(across(taget_col, as.numeric))
pheno_data_filter$BMI <- ifelse(pheno_data_filter$BMI >= 25, "High", "Low")

# Summarize the phenotype data
library(gtsummary)
tbl_summary(data = pheno_data_filter |>
              dplyr::select(all_of(var_focus_3)),
            by = "BMI",
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{median} ({min}, {max})", "{mean}+/-{sd}")),
            missing_text = "Missing, n") |>
  modify_table_body(
    dplyr::mutate,
    label = case_when(label == "Median (Range)" ~ "Median (range)",
                      label == "Mean+/-SD" ~ "Mean+/-sd",
                      TRUE ~ label)) |>
  modify_footnote(update = everything() ~ NA) |>
  bold_labels() |>
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3),
        test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 5000)) 

