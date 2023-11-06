# **Distribution**
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

# Calculate the probability of obtaining a normalized gene expression value that is within three standard deviations of the mean for gene_expression ∼ N(0, 1)
pnorm(q = 3, mean = 0, sd = 1) - pnorm(q = -3, mean = 0, sd = 1)
# [1] 0.9973002

# What is the probability that out of a sample size of 1,000 individuals who have taken a new experimental drug, we observe 0.38 or more experiencing significant side effects? 
# Use a normal approximation to the binomial to estimate this probability and compare it to the exact probability from a binomial distribution.
1 - pnorm(0.379, mean = pi_true, sd = sqrt(pi_true * (1 - pi_true) / 1000))
1 - pbinom(0.379 * 1000, size = 1000, prob = pi_true)

# Estimate the kurtosis
library(haven)
library(e1071)
acath_data <- read_dta("acath.dta")
kurtosis(acath_data$cad_dur)

# Test if the sample mean is equal to 280 (population mean)
library(BSDA)
z.test(x = rep(275, 50), 
       alternative = "two.sided", 
       mu = 280, 
       sigma.x = 10, 
       conf.level = 0.95)
p_val = 2*pnorm(275, 280, 10/sqrt(50))

# Calculate the log-likelihood of a binomial parameter
# Method 1:
n <- 100
x <- 30
pi <- seq(0, 1, by = 0.01)
log_likelihood <- log(choose(n, x)) + x * log(pi) + (n - x) * log(1 - pi)
# Method 2:
dbinom(30, 100, prob = seq(0, 1, by = 0.01), log = TRUE) # log = F -> estimate the likelihood

# Product of the vector elements
prod(1:7)
# [1] 5040

# **Matrix**
# Calculate the Moore-Penrose generalized inverse of a matrix X
ginv(t(x) %*% x)

# **Model**
# Fit the generalized models
mdl.fit <- pscl::zeroinfl(count~1|1, link = "logit", dist = "negbin", data = test.dat)
mdl.fit <- VGAM::vglm(count~1, negbinomial, data = test.dat)

# Build the Lasso and Elastic-Net regularized generalized linear models
# https://glmnet.stanford.edu/articles/glmnet.html
