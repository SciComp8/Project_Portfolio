library(ggplot2)
organized_duplicate_rate_data[,variable:=factor(variable,levels=c("BMI", "AGE","SEX", "MHABNWBC", "BMIxSEX"))]
organized_duplicate_rate_data[,cDEGs.duplicate.rate.median:=as.numeric(cDEGs.duplicate.rate.median)]
organized_duplicate_rate_data[,cDEGs.duplicate.rate.mean:=as.numeric(cDEGs.duplicate.rate.mean)]
date.analysis <- format(Sys.Date(), "%Y%b%d")
g1 <- organized_duplicate_rate_data |>
  ggplot(mapping = aes(y = cDEGs.duplicate.rate.median, x = method, group = method)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = method)) + 
  facet_wrap(vars(variable)) + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     breaks = seq(0, 1, 0.1)) + 
  labs(y = "Median duplicate rate of cDEGs", x = "Method") + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5))
ggsave(filename = sprintf("../ApplicationResult/Multi_Interaction/RandomSeed/DuplicatedRateMatrix/Boxplot/%s_median.png", date.analysis),
       plot = g1, device = "png", width = 9, height = 6, units = "in")

g2 <- organized_duplicate_rate_data |>
  ggplot(mapping = aes(y = cDEGs.duplicate.rate.mean, x = method, group = method)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = method)) + 
  facet_wrap(vars(variable)) + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     breaks = seq(0, 1, 0.1)) + 
  labs(y = "Mean duplicate rate of cDEGs", x = "Method") + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5))
ggsave(filename = sprintf("../ApplicationResult/Multi_Interaction/RandomSeed/DuplicatedRateMatrix/Boxplot/%s_mean.png", date.analysis),
       plot = g2, device = "png", width = 9, height = 6, units = "in")
