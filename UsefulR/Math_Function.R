# Calculate the Moore-Penrose generalized inverse of a matrix X
ginv(t(x) %*% x)

# Calculate the distribution function for the F distribution with df1 and df2 degrees of freedom
pf(F.stat, df.M, df.E, lower.tail = F)

# Calculate and visualize the probability density function and cumulative distribution function for the normal distribution
x <- seq(0, 40, by = 0.1)
mean_value <- 25
sd_value <- 5
pdf_values <- dnorm(x, mean = mean_value, sd = sd_value)
plot(x, pdf_values, type = "l", col = "blue", lwd = 2, xlab = "X-axis", ylab = "PDF", main = "PDF of Normal Distribution")
# type
# "p" for points,
# "l" for lines,
# "b" for both,
# "c" for the lines part alone of "b",
# "o" for both ‘overplotted’,
# "h" for ‘histogram’ like (or ‘high-density’) vertical lines,
# "s" for stair steps,
# "S" for other steps, see ‘Details’ below,
# "n" for no plotting.
grid() # Add gridlines (optional)

cdf_values <- pnorm(x, mean = mean_value, sd = sd_value)
plot(x, pdf_values, type = "l", col = "orange", lwd = 2, xlab = "X-axis", ylab = "PDF", main = "CDF of Normal Distribution")
grid()

# Fit the generalized models
mdl.fit <- pscl::zeroinfl(count~1|1, link = "logit", dist = "negbin", data = test.dat)
mdl.fit <- VGAM::vglm(count~1, negbinomial, data = test.dat)
