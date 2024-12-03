# References:
# https://github.com/SciComp8/Project_Programming/blob/main/UsefulR/attach.R
# Make a density plot
plot(density(df$var), xlab="RNA velocity", main="Distribution of RNA velocity",
     lwd=3, col="chocolate", ylim=c(0, 0.08))

dat0 <- read_xls(
  path = "../Data/Raw/41586_2012_BFnature11735_MOESM103_ESM.xls", 
  sheet = "Table S11",
  col_names = T,
  na = c("Unknown", "unknown", "N/A", "n/a"),
  skip = 1
) 
dat0 <- dat0[, -13]
dat0 <- dat0 |>
  clean_names() 
dim(dat0) 

summary(dat0$act_area_raw)
plot(density(dat0$act_area_raw), xlab="Raw Activity Area", main="Distribution of Raw Activity Area",
     lwd=3, col="chocolate")

# Make a histogram
hist(df$var, breaks=50, xlab="RNA velocity", main="Distribution of RNA velocity")
hist(peak_data$peak_ratio, breaks=50, xlab="Peak ratio (contaminated peak / normal peak)", main="Distribution of peak ratio")
hist_data <- hist(peak_data$peak_ratio, plot = FALSE) # Generate histogram data without plotting
hist_data$counts <- log2(hist_data$counts+1) # Transform frequencies to log2 scale; +1 to avoid log(0); 
# Plot the histogram with log2-transformed frequencies
plot(hist_data, freq = TRUE, breaks=50, 
     main = "Distribution of peak ratio with log2 transformed frequencies",
     ylab = "log2 (rrequency)",
     xlab = "Peak ratio",
     ylim = c(0, max(hist_data$counts))) # Adjust ylim based on transformed frequencies

# Quick code using base R functions
airquality #dataset from New York Air Quality Measurements
attach(airquality) 
hist(Temp)
pairs(airquality[, 1:4]) #plot all possible pairwise scatterplots
