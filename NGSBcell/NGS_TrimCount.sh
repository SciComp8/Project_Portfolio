#!/bin/bash

# Set path variables
dir=/path/to/working/directory 
fastq_dir=$dir/data/samples
trim_dir=$dir/trimmed_reads
kallisto_dir=$dir/kallisto_results
DESeq2_dir=$dir/DESeq2_results
transcriptome_idx=$dir/mus_musculus/transcriptome.idx

# Retrieve sample names, conditions, and number of biological replicates
sample_name=( $(ls $fastq_dir/*_R1_001_copy1.fastq.gz | awk -F/ '{print $NF}' | awk -F_ '{print $1}' | sort -u) )
condition=( $(echo ${sample_name[@]} | tr ' ' '\n' | awk -F_ '{print $1}' | sort -u) )
bio_replicate=( $(echo ${condition[@]} | tr ' ' '\n' | uniq -c | awk '{print $1}') )

# Construct parameter strings for Snakemake rule
sample_list=$(printf '"%s",' "${sample_name[@]}")
sample_list=${sample_list%,}
condition_list=$(printf '"%s",' "${condition[@]}")
condition_list=${condition_list%,}
condition_times=$(printf '%s,' "${bio_replicate[@]}")
condition_times=${condition_times%,}

# Run Snakemake rule
snakemake --config sample_list=$sample_list condition_list=$condition_list condition_times=$condition_times

# Trim reads using Trim Galore
for sample in ${sample_name[@]}; do
  for replicate in {1..${bio_replicate[0]}}; do
    input=$fastq_dir/${sample}_R1_001_copy${replicate}.fastq.gz
    output=$trim_dir/$sample/${sample}_R1_001_copy${replicate}_trimmed.fq.gz
    trim_galore $input --cores 4 --quality 20 --output_dir $trim_dir/$sample --no_report_file
  done
done

# Quantify gene expression using Kallisto
for sample in ${sample_name[@]}; do
  input=$(ls $trim_dir/$sample/*_trimmed.fq.gz | sort)
  output_dir=$kallisto_dir/$sample
  kallisto quant --index $transcriptome_idx --plaintext --output-dir $output_dir --single -l 200 -s 20 $input
done

# Generate count matrix
cd $kallisto_dir
mkdir -p count_matrix
for sample in ${sample_name[@]}; do
  mkdir -p count_matrix/$sample
  kallisto-to-matrix -o count_matrix/$sample abundance.tsv
done

# Perform differential gene expression analysis using DESeq2
cd $DESeq2_dir
Rscript -e $sample_list $condition_list $condition_times 'library(rmarkdown); library(rmdformats); rmarkdown::render(input = "NGS_DEReport.Rmd", output_format = "all", output_file = "Final Report")'
