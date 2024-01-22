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

