vim fastq2bam.slurm

# Inside fastq2bam.slurm
#!/bin/sh
#SBATCH --partition=panda
#SBATCH --time=00:20:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=10

spack load samtools@1.13%gcc@8.2.0
spack load /ztcq4ql

bowtie2 -p 10 -q --local -x /home/your_name/project/index/hg38_chr1 -U /home/your_name/project/SRRXXXX.fastq -S /home/your_name/project/SRRXXXX.sam 2> /home/your_name/project/SRRXXXX.log | samtools view -@ 10 -Sb -o /home/your_name/project/SRRXXXX.bam 1> /home/your_name/project/SRRXXXX_summary.txt 2> /home/your_name/project/SRRXXXX_error.txt & 

# Outside fastq2bam.slurm
sbatch fastq2bam.slurm
