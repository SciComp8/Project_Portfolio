# Calculate the Moore-Penrose generalized inverse of a matrix X
ginv(t(x) %*% x)

# Calcuate the distribution function for the F distribution with df1 and df2 degrees of freedom
pf(Fstatistic, df_M, df_E, lower.tail = F)

# Fit the generalized models
mdl.fit <- pscl::zeroinfl(count~1|1, link="logit", dist="negbin",data=test.dat)
mdl.fit <- VGAM::vglm(count~1, negbinomial, data=test.dat)
