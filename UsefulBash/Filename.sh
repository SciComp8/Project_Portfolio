# Extract the base filenames (without the .fastq.gz extension) of the input files specified by the first and second arguments 
r1=$(basename $1 .fastq.gz)
r2=$(basename $2 .fastq.gz)
