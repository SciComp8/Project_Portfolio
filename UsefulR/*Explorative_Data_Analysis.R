# Make a density plot
plot(density(df$var), xlab="RNA velocity", main="Distribution of RNA velocity",
     lwd=3, col="chocolate", ylim=c(0, 0.08))

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
