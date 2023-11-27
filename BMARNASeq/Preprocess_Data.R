#----Load data files----
library(readxl)
library(tidyverse)
library(gtsummary)
# BiocManager::install(version = "3.18")
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

data_dir <- '../data/' # >
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

#----Clean data----
names(pheno_data)
table(pheno_data$SMTSD, useNA='ifany') # View the distribution of tissue type
idx_sub_fat <- which(pheno_data$SMTSD=='Adipose - Subcutaneous')
dict_name <- dict_full$name
length(dict_name); ncol(pheno_data)
colnames(pheno_data)[!colnames(pheno_data) %in% dict_name]
# [1] "dbGaP_Subject_ID.x" "dbGaP_Subject_ID.y" "dbGaP_Sample_ID" 
pheno_data_filter <- pheno_data[, colnames(pheno_data) %in% dict_name]
var_str <- dict_full$type=='string'
sum(var_str) # Number of string variables
var_focus <- c('COHORT', as.character(dict_full$name[!var_str]))
var_cat <- c(1, as.numeric(dict_full$type[!var_str] == 'integer, encoded value'))
colSums(is.na(pheno_data_filter[, var_focus]))
var_remove <- which(colSums(is.na(pheno_data_filter[, var_focus])) > 70) # ? 70
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
  mutate(across(all_of(taget_col), as.numeric))
pheno_data_filter$BMI <- ifelse(pheno_data_filter$BMI >= 25, "High", "Low")

#----Take a view at the distribution of the phenotypes---
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

# >
pheno_data_final <- pheno_data_filter[, var_focus_3]  
dim(pheno_data_final)
rownames(pheno_data_final) <- pheno_data$SAMPID
rownames(pheno_data_final) <- gsub(pattern = "-", replacement = ".", x = rownames(pheno_data_final))
expr_data <- read.csv(paste0(data_dir, 'expr_data.csv'))
ensg_id <- rownames(expr_data)
dim(expr_data)
if (identical(rownames(pheno_data_final), colnames(expr_data))) {
  print('The elements in one vector are in the same order as they appear in another vector!')
  expr_data_filter <- expr_data
} else {
  idx_sample <- match(rownames(pheno_data_final), colnames(expr_data))
  # match returns a vector of the positions of (first) matches of its first argument in its second.
  expr_data_filter <- expr_data[, idx_sample]
}

# The following 2 code lines are equivalent
# identical(rownames(pheno_data_final), colnames(expr_data_filter))
# all(match(rownames(pheno_data_final), colnames(expr_data)) == 1:nrow(pheno_data_final))

colnames(expr_data_filter) # These are sample IDs
dim(expr_data_filter) 

# Build the DESeqDataSet
pheno_data_final$BMI <- as.factor(pheno_data_final$BMI)
dds <- DESeqDataSetFromMatrix(countData = expr_data_filter,
                              colData = pheno_data_final,
                              design = ~ BMI)
dds
dds <- DESeq(dds)
dds
vsd <- vst(dds, blind=FALSE)

#----Evaluate the data quality by sample clustering and visualization----
sample_distance <- dist(t(assay(vsd)), method='euclidean') 
# dist: This function computes and returns the distances between the rows of a data matrix.
sample_distance_mat <- as.matrix(sample_distance)
row_order <- order(paste0(vsd$BMI, '-', rownames(sample_distance_mat)), decreasing=T)
column_order <- order(paste0(vsd$BMI, '-', colnames(sample_distance_mat)), decreasing=T)
new_rownames <- rownames(sample_distance_mat)[row_order]
new_columnnames <- colnames(sample_distance_mat)[column_order]
sample_distance_mat_2 <- sample_distance_mat[match(new_rownames, rownames(sample_distance_mat)), match(new_columnnames, colnames(sample_distance_mat))] 

colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)
# pheatmap(mat=sample_distance_mat,
#          clustering_distance_rows=sample_distance,
#          clustering_distance_cols=sample_distance,
#          col=colors)

# >
num_blocks <- nrow(sample_distance_mat) %/% 6
for (i in 1:num_blocks) {
  print(paste0('Plotting ', i, ' heatmap'))
  pheatmap(mat=sample_distance_mat_2[((i-1)*6+1):(i*6), ((i-1)*6+1):(i*6)],
           col=colors)
}

# >
plotPCA(vsd, intgroup=c('BMI')) +
  labs(color='BMI')

# > TD
select_top30 <- order(rowMeans(counts(dds,normalized=T)), decreasing=T)[1:30]
pheatmap(assay(vsd)[select_top30,], 
         cluster_rows=F, show_rownames=F, cluster_cols=F, 
         annotation_col=as.data.frame(pheno_data_filter[,'BMI']),
         annotation_colors = list(BMI = c("Low" = "blue", "High" = "red")))

# > 
data_dir <- '../data/'
pheno_data <- readRDS(paste0(data_dir, 'pheno_data.RDS'))
expr_data <- readRDS(paste0(data_dir, 'expr_data.RDS'))

#----Remove missing values in the phenotype data---- 
data_dir <- '../data/'
pheno_data <- readRDS(paste0(data_dir, 'pheno_data.RDS'))
expr_data <- readRDS(paste0(data_dir, 'expr_data.RDS'))

var_13 <- c('SEX', 'AGE', 'BMI', 'MHABNWBC', 'MHARTHTS', 
            'MHASTHMA', 'MHBCTINF', 'MHBLDDND', 'MHCOCAINE5', 
            'MHCOPD', 'MHCVD', 'MHDRNKSTS', 'SMTSPAX')

# MHABNWBC: normal or abnormal white blood cell
# MHARTHTS: has arthritis or not
# MHASTHMA: has asthma or not 
# MHBCTINF: has bacterial infections or not
# MHBLDDND: has blood donation been denied in the past or not
# MHCOCAINE5: use cocaine in the past 5 years or not
# MHCOPD: has chronic respiratory disease or not
# MHCVD: has cerebrovascular disease disorder or not
# MHDRNKSTS: drink or not
# SMTSPAX: time a sample spent in the PAXgene fixative

pheno_data_13 <- pheno_data[, var_13]
all_char_var <- sapply(names(pheno_data_13)[sapply(pheno_data_13, is.character)], function(x) with(pheno_data_13, table(get(x), useNA = "ifany"))) 
all_num_var <- sapply(names(pheno_data_13)[sapply(pheno_data_13, is.numeric)], function(x) with(pheno_data_13, table(get(x), useNA = "ifany"))) 
idx_exclude <- c(which(pheno_data_13$MHABNWBC==99), 
                 which(pheno_data_13$MHBCTINF==99), 
                 which(pheno_data_13$MHBLDDND==99), 
                 which(pheno_data_13$MHCOCAINE5==99), 
                 which(pheno_data_13$MHCOPD==99),
                 which(pheno_data_13$MHCVD==99),
                 which(is.na(pheno_data_13$MHDRNKSTS))) |> unique()
length(idx_exclude) # 31

pheno_data_13 <- pheno_data_13[-idx_exclude, ]  
expr_data_13 <- expr_data[, -idx_exclude]    
dim(pheno_data_13) 
dim(expr_data_13) 
sum(is.na(pheno_data_13)) # 0

#----Exclude genes from the RNA-seq data which have Counts per Million greater than 0.2 but appear in 10 or fewer samples---- 
num_0.2 <- rowSums(cpm(expr_data_13) > 0.2)
# expr_cpm <- cpm(expr_data_13); num_0.2 <- apply(expr_cpm, 1, function(x) sum(x > 0.2))
expr_data_13 <- expr_data_13[num_0.2 > 10, ]
length(which(num_0.2 <= 10)) 
dim(expr_data_13) 

#----Recode values in the phenotype data----
for (col in var_13) {
  if (col %in% c("SEX", "AGE", "SMTSPAX", "MHDRNKSTS")) {
    # Handle special cases
    if (col == "SEX") {
      pheno_data_13[[col]] <- factor(ifelse(pheno_data_13[[col]] == 1, "male", "female"), levels = c("female", "male"))
    } else if (col == "AGE") {
      age_median <- median(pheno_data_13[[col]])
      pheno_data_13[[col]] <- factor(ifelse(pheno_data_13[[col]] < age_median, "young", "old"), levels = c("young", "old"))
    } else if (col == "SMTSPAX") {
      smtspax_median <- median(pheno_data_13[[col]])
      pheno_data_13[[col]] <- factor(ifelse(pheno_data_13[[col]] < smtspax_median, "low", "high"), levels = c("low", "high"))
    } else if (col == "MHDRNKSTS") {
      pheno_data_13[[col]] <- factor(ifelse(pheno_data_13[[col]] == "No", "no", "yes"), levels = c("no", "yes"))
    }
  } else {
    # Handle general cases
    pheno_data_13[[col]] <- factor(ifelse(pheno_data_13[[col]] == 0, "no", "yes"), levels = c("no", "yes"))
  }
}

all_factor_var <- sapply(names(pheno_data_13)[sapply(pheno_data_13, is.factor)], function(x) with(pheno_data_13, table(get(x), useNA = "ifany")), simplify = F) 
all_num_var <- sapply(names(pheno_data_13)[sapply(pheno_data_13, is.numeric)], function(x) with(pheno_data_13, table(get(x), useNA = "ifany"))) 

saveRDS(pheno_data_13, file=paste0(data_dir, 'dat.pheno.Subcutaneous'))
saveRDS(expr_data_13, file=paste0(data_dir, 'dat.expr.Subcutaneous'))
