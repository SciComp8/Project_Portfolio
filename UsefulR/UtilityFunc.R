# Deletes the file(s) or directories 
unlink("~/Downloads/Lab7_GWAS/lab7_cache", recursive = T)

# Read the file
read.delim(file, sep = ",")

# Stack two columns into one long vector == translocate the second column into the bottom of the first column
mapply(c, geno_import[, seq(1, ncol(geno_import), 2)], geno_import[, seq(2, ncol(geno_import), 2)])

# Check the arguments of a function
args(mapply)

# Calculate the Moore-Penrose generalized inverse of a matrix X
ginv(t(x) %*% x)

# Calcuate the distribution function for the F distribution with df1 and df2 degrees of freedom
pf(Fstatistic, df_M, df_E, lower.tail = F)
