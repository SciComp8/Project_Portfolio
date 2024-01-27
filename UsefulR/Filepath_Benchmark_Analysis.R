out_file_info <- function(file_info_table) {
  unique_conditions <- unique(file_info_table$conditions)
  condition_1 <- unique_conditions[1]
  condition_2 <- unique_conditions[2]
  
  samples_1 <- file_info_table$samples[file_info_table$conditions == condition_1]
  samples_2 <- file_info_table$samples[file_info_table$conditions == condition_2]
  
  extract_dir <- function(condition, column_name, suffix) {
    files <- file_info_table[file_info_table$conditions == condition, column_name]
    if (column_name %in% c("files_bam", "dirs_salmon")) {
      setNames(files, get(sprintf("samples_%s", suffix)))
    } else {
      setNames(files, get(sprintf("samples_%s", suffix))) |> sapply(function(x) unlist(strsplit(x, ";"))[suffix])
    }
  }
  
  bams_1 <- extract_dir(condition_1, "files_bam", 1)
  bams_2 <- extract_dir(condition_2, "files_bam", 2)
  
  if ("dirs_salmon" %in% colnames(file_info_table)) {
    salmons_1 <- extract_dir(condition_1, "dirs_salmon", 1)
    salmons_2 <- extract_dir(condition_2, "dirs_salmon", 2)
  } else {
    salmons_1 <- salmons_2 <- NULL
  }
  
  if ("files_fq" %in% colnames(file_info_table)) {
    fqs_R1_1 <- extract_dir(condition_1, "files_fq", 1)
    fqs_R2_1 <- extract_dir(condition_1, "files_fq", 1)
    fqs_R1_2 <- extract_dir(condition_2, "files_fq", 2)
    fqs_R2_2 <- extract_dir(condition_2, "files_fq", 2)
  } else {
    fqs_R1_1 <- fqs_R2_1 <- fqs_R1_2 <- fqs_R2_2 <- NULL
  }
  
  info_list <- list(
    condition_1 = condition_1,
    condition_2 = condition_2,
    samples_1 = samples_1,
    samples_2 = samples_2,
    bams_1 = bams_1,
    bams_2 = bams_2,
    salmons_1 = salmons_1,
    salmons_2 = salmons_2,
    fqs_R1_1 = fqs_R1_1,
    fqs_R2_1 = fqs_R2_1,
    fqs_R1_2 = fqs_R1_2,
    fqs_R2_2 = fqs_R2_2
  )
  
  return(info_list)
}

### FASTQ files ###
sample_group_1 <- paste0("sample", 1:4)
sample_group_2 <- paste0("sample", 5:8)
suffix_fq <- "fq.gz"
dir_fq <- "~/projects/bulkRNAseq/splicing/data/fq"

(fq_group_1 <- paste(
  sprintf("%s/%s_R1.%s", dir_fq, sample_group_1, suffix_fq),
  sprintf("%s/%s_R2.%s", dir_fq, sample_group_1, suffix_fq),
  sep = ";"
))
# [1] "~/projects/bulkRNAseq/splicing/data/fq/sample1_R1.fq.gz;~/projects/bulkRNAseq/splicing/data/fq/sample1_R2.fq.gz"
# [2] "~/projects/bulkRNAseq/splicing/data/fq/sample2_R1.fq.gz;~/projects/bulkRNAseq/splicing/data/fq/sample2_R2.fq.gz"
# [3] "~/projects/bulkRNAseq/splicing/data/fq/sample3_R1.fq.gz;~/projects/bulkRNAseq/splicing/data/fq/sample3_R2.fq.gz"
# [4] "~/projects/bulkRNAseq/splicing/data/fq/sample4_R1.fq.gz;~/projects/bulkRNAseq/splicing/data/fq/sample4_R2.fq.gz"

fq_group_2 <- paste(
  sprintf("%s/%s_R1.%s", dir_fq, sample_group_2, suffix_fq),
  sprintf("%s/%s_R2.%s", dir_fq, sample_group_2, suffix_fq),
  sep = ";"
)

### BAM files ###
suffix_bam <- "SortedByCoord.bam"
dir_bam <- "~/projects/bulkRNAseq/splicing/data/star"

(bam_group_1 <- sprintf("%s/%s/%s.%s", dir_bam, sample_group_1, sample_group_1, suffix_bam))
# [1] "~/projects/bulkRNAseq/splicing/data/star/sample1/sample1.SortedByCoord.bam"
# [2] "~/projects/bulkRNAseq/splicing/data/star/sample2/sample2.SortedByCoord.bam"
# [3] "~/projects/bulkRNAseq/splicing/data/star/sample3/sample3.SortedByCoord.bam"
# [4] "~/projects/bulkRNAseq/splicing/data/star/sample4/sample4.SortedByCoord.bam"

bam_group_2 <- sprintf("%s/%s/%s.%s", dir_bam, sample_group_2, sample_group_2, suffix_bam)

### Generate a summary table of input files ###
file_info_table <- data.frame(
  samples = c(sample_group_1, sample_group_1),
  conditions = c(rep("wt", 4), rep("kd", 4)),
  files_bam = c(bam_group_1, bam_group_2),
  files_fq = c(fq_group_1, fq_group_2)
)  

### Output the meta information of input files ###
out_file_info(file_info_table)

# $condition_1
# [1] "wt"
# 
# $condition_2
# [1] "kd"
# 
# $samples_1
# [1] "sample1" "sample2"
# [3] "sample3" "sample4"
# 
# $samples_2
# [1] "sample1" "sample2"
# [3] "sample3" "sample4"
# 
# $bams_1
# sample1 
# "~/projects/bulkRNAseq/splicing/data/star/sample1/sample1.SortedByCoord.bam" 
# sample2 
# "~/projects/bulkRNAseq/splicing/data/star/sample2/sample2.SortedByCoord.bam" 
# sample3 
# "~/projects/bulkRNAseq/splicing/data/star/sample3/sample3.SortedByCoord.bam" 
# sample4 
# "~/projects/bulkRNAseq/splicing/data/star/sample4/sample4.SortedByCoord.bam" 
# 
# $bams_2
# sample1 
# "~/projects/bulkRNAseq/splicing/data/star/sample5/sample5.SortedByCoord.bam" 
# sample2 
# "~/projects/bulkRNAseq/splicing/data/star/sample6/sample6.SortedByCoord.bam" 
# sample3 
# "~/projects/bulkRNAseq/splicing/data/star/sample7/sample7.SortedByCoord.bam" 
# sample4 
# "~/projects/bulkRNAseq/splicing/data/star/sample8/sample8.SortedByCoord.bam" 
# 
# $salmons_1
# NULL
# 
# $salmons_2
# NULL
# 
# $fqs_R1_1
# sample1 
# "~/projects/bulkRNAseq/splicing/data/fq/sample1_R1.fq.gz" 
# sample2 
# "~/projects/bulkRNAseq/splicing/data/fq/sample2_R1.fq.gz" 
# sample3 
# "~/projects/bulkRNAseq/splicing/data/fq/sample3_R1.fq.gz" 
# sample4 
# "~/projects/bulkRNAseq/splicing/data/fq/sample4_R1.fq.gz" 
# 
# $fqs_R2_1
# sample1 
# "~/projects/bulkRNAseq/splicing/data/fq/sample1_R1.fq.gz" 
# sample2 
# "~/projects/bulkRNAseq/splicing/data/fq/sample2_R1.fq.gz" 
# sample3 
# "~/projects/bulkRNAseq/splicing/data/fq/sample3_R1.fq.gz" 
# sample4 
# "~/projects/bulkRNAseq/splicing/data/fq/sample4_R1.fq.gz" 
# 
# $fqs_R1_2
# sample1 
# "~/projects/bulkRNAseq/splicing/data/fq/sample5_R2.fq.gz" 
# sample2 
# "~/projects/bulkRNAseq/splicing/data/fq/sample6_R2.fq.gz" 
# sample3 
# "~/projects/bulkRNAseq/splicing/data/fq/sample7_R2.fq.gz" 
# sample4 
# "~/projects/bulkRNAseq/splicing/data/fq/sample8_R2.fq.gz" 
# 
# $fqs_R2_2
# sample1 
# "~/projects/bulkRNAseq/splicing/data/fq/sample5_R2.fq.gz" 
# sample2 
# "~/projects/bulkRNAseq/splicing/data/fq/sample6_R2.fq.gz" 
# sample3 
# "~/projects/bulkRNAseq/splicing/data/fq/sample7_R2.fq.gz" 
# sample4 
# "~/projects/bulkRNAseq/splicing/data/fq/sample8_R2.fq.gz" 
