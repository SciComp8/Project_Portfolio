library(mixOmics)
dat <- linnerud[["physiological"]]
cov(dat)
cor(dat, method='pearson')
