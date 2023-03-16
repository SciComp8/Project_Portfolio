# Deletes the file(s) or directories 
unlink("~/Downloads/Lab7_GWAS/lab7_cache", recursive = T)

# Read the file
read.delim(file, sep = ",")


# Stack two columns into one long vector == translocate the second column into the bottom of the first column
mapply(c, geno_import[, seq(1, ncol(geno_import), 2)], geno_import[, seq(2, ncol(geno_import), 2)])
