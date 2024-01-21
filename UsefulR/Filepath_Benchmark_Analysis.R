# Define variables for sample names and suffix
sample_group_1 <- paste0("sample", 1:4)
sample_group_2 <- paste0("sample", 5:8)
suffix_bam <- "SortedByCoord.bam"
dir_bam <- "~/projects/bulkRNAseq/splicing/data/star/"

# Create file paths for group 1 samples
bam_group_1 <- sprintf("%s/%s/%s.%s", dir_bam, sample_group_1, sample_group_1, suffix_bam)

# Create file paths for group 2 samples
bam_group_2 <- sprintf("%s/%s/%s.%s", dir_bam, sample_group_2, sample_group_2, suffix_bam)
