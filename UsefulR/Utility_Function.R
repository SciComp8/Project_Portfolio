# Deletes the file(s) or directories 
unlink("~/Downloads/Lab7_scRNAseq/lab7_cache", recursive = T)

# Construct the path
file.path("","p1","p2","p3", "filename")
[1] "/p1/p2/p3/filename"

# Extract the components from a path
dirname(file.path("","p1","p2","p3", "filename"))
[1] "/p1/p2/p3"
basename(file.path("","p1","p2","p3", "filename"))
[1] "filename"

# Download the file
download.file(url = "https://XXX.tsv", destfile = "XXX.csv")

# Read the file
read.delim(file, sep = ",")
readr::read_tsv("XXX.tsv")

# Write the file
file.name = "Matrix.txt"
write(sprintf("probability matrix: x= %s", ncol(object1)),
      file.name, append = T)
write.table(object2,
            file.name, sep = "\t",
            quote = F, 
            col.names = F,
            append = T)

# Stack two columns into one long vector == translocate the second column into the bottom of the first column
mapply(c, geno_import[, seq(1, ncol(geno_import), 2)], geno_import[, seq(2, ncol(geno_import), 2)])

# Check the arguments of a function
args(mapply)
