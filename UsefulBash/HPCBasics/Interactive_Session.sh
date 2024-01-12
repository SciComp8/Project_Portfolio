# Reference: 
# https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome
# https://stackoverflow.com/questions/38502194/could-not-locate-a-bowtie-index-corresponding-to-basename
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

# In local terminal
mv Downloads/chr1.fa.gz ~
gzip -dk chr1.fa.gz
rsync -avP /Users/your_name/SRRXXX.fastq /Users/your_name/chr1.fa username@aphrodite.med.cornell.edu:/home/username/project/index

# In HPC
mkdir -p project/index
spack load /ztcq4ql
bowtie2-build index/chr1.fa hg38_chr1 --threads 10

bowtie2 -p 10 -q --local -x /home/username/project/index/hg38_chr1 -U /home/username/project/SRRXXX.fastq -S /home/username/project/SRRXXX.sam 2> ./SRRXXX.log
cat SRRXXX.log
