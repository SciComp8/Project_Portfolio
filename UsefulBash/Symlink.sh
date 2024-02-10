mkdir Documents/fastq && mkdir Documents/data && ln -s Documents/fastq Documents/data
find ~/Documents -type l
ls -l
ln -fs Documents/fastq/* Documents/data
