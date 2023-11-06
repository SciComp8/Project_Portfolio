# Quick way
qqnorm(data_vec, main="Normal Q-Q plot for the sample variance of sampling distribution", col="midnightblue")
qqline(data_vec, col="indianred", lwd=2)

# QQ plot shared with collaborators
ggplot() +
  stat_qq(aes(sample=data_vec)) +
  stat_qq_line(aes(sample=data_vec)) +
  labs(x="Theoretical Quantiles", y="Sample Quantiles", title="QQ Plot") +
  theme_bw()
