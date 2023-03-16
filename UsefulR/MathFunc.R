# Calculate the Moore-Penrose generalized inverse of a matrix X
ginv(t(x) %*% x)

# Calcuate the distribution function for the F distribution with df1 and df2 degrees of freedom
pf(Fstatistic, df_M, df_E, lower.tail = F)
