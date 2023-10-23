# Batch view the filenames
for filename in relative_path/*.fastq.gz; do echo $filename; done

# Extract the base filenames (without the .tsv.gz extension) of the input files specified by the first and second arguments 
r1=$(basename $1 .tsv.gz)
r2=$(basename $2 .tsv.gz)
echo $r1 $r2
# bash Filename.sh barcodes1.tsv.gz barcodes2.tsv.gz


# Remove files with specific shared base filenames
# More info: https://unix.stackexchange.com/questions/306940/what-is-the-purpose-of-the-do-keyword-in-bash-for-loops
for i in `ls | grep .fastq`
do
  basename=$(echo $i | sed 's/\50K_2_bowtie2//g')
  rm $basename.sam
done

# Extract the specific fields in the string text and concatenate the output with no separator in between
filename="SRR8990876_1.paired.fastq.gz"
echo $filename | awk -F"." 'BEGIN{OFS="_"} {print $1,$2}'
# Output: SRR8990876_1_paired
echo $filename | awk -F"." -v OFS="_" '{print $1,$2}'
# Output: SRR8990876_1_paired


# Output variables into a file with the specific filename
filename=$(echo $(basename barcodes1 .tsv.gz))
for i in `echo $var1`
do
  echo $var1 $var2 $var3 $i >>$filename.txt 
done


# Compress/decompress a file with the specific filename
filename=$(echo $(basename SRR8990876_1.paired .fastq.gz))
gzip $filename.fastq
gunzip $filename.fastq.gz


# Create a backup copy of the file with the specified filename
fastq=$1 # Store the first argument passed to the program as variable fastq
cp $fastq rna_copy.fastq

