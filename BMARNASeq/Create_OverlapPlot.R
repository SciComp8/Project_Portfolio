load("../ApplicationData/derived/overlapplot_image.RData")

suppressPackageStartupMessages(easypackages::libraries("microbenchmark", "bannerCommenter", "fst", "parallel", "ggplot2", "data.table"))

## Batch plotting and saving - num.cDEGs
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.vec <- c("BMI", "AGE", "SEX", "WBC", "BMIxSEX")

for (var.i in var.vec) {
  g <- boxplot.data.new[variable==var.i] |>
    ggplot(aes(y = num.cDEGs, x = method, group = method, fill = method)) + 
    geom_boxplot() + 
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "#0073C2FF") +  
    facet_wrap(vars(threshold)) + 
    labs(y = "Number of common differentially expressed genes", x = "Method") + 
    theme_bw() + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5))
    ggsave(filename = sprintf("../ApplicationResult/Multi_Interaction/RandomSeed/numcDEGs/Boxplot/%s_%s.png", date.analysis, var.i),
           plot = g, device = "png", width = 9, height = 6, units = "in")
}

## Batch plotting and saving - rate.cDEGs
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.vec <- c("BMI", "AGE", "SEX", "WBC", "BMIxSEX")

for (var.i in var.vec) {
  g <- boxplot.data.new2[variable==var.i] |>
    ggplot(aes(y = rate.cDEGs, x = method, group = method, fill = method)) + 
    geom_boxplot() + 
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "#0073C2FF") +  
    facet_wrap(vars(threshold)) + 
    labs(y = "Rate of common differentially expressed genes", x = "Method") + 
    theme_bw() + 
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5))
    ggsave(filename = sprintf("../ApplicationResult/Multi_Interaction/RandomSeed/ratecDEGs/Boxplot/%s_%s.png", date.analysis, var.i),
           plot = g, device = "png", width = 9, height = 6, units = "in")
}
