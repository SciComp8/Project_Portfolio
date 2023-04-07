# Extract the base filenames (without the .tsv.gz extension) of the input files specified by the first and second arguments 
r1=$(basename $1 .tsv.gz)
r2=$(basename $2 .tsv.gz)
echo $r1 $r2
# bash Filename.sh barcodes1.tsv.gz barcodes2.tsv.gz
