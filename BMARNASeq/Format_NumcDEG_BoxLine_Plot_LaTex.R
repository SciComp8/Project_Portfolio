source("theme_BMA.R")
library(cowplot)

##------Preprocess data------
cDEGs_num_uni <- readRDS("../ApplicationResult/Uni/RandomSeed/numcDEGs/Boxplot/boxplot.data.num.BMA.RDS")
cDEGs_num_uni <- cDEGs_num_uni[cDEGs_num_uni$method != "BMAseq", ]
cDEGs_num_multi <- readRDS("../ApplicationResult/Multi/RandomSeed/numcDEGs/Boxplot/boxplot.data.num.RDS")

cDEGs_num_uni$method <- paste0(cDEGs_num_uni$method, "_UVM")
cDEGs_num_multi$method <- ifelse(cDEGs_num_multi$method != "BMAseq", paste0(cDEGs_num_multi$method, "_MVM"), cDEGs_num_multi$method)
cDEGs_num_all <- rbind(cDEGs_num_uni, cDEGs_num_multi)
cDEGs_num_all$method <- factor(cDEGs_num_all$method,
                               levels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM"))
cDEGs_num_all$variable <- factor(cDEGs_num_all$variable, 
                                 levels = c("BMI", "AGE", "SEX", "MHABNWBC"))
cDEGs_num_all$random.seed <- factor(cDEGs_rate_all$random.seed,
                                    levels = c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999))

date.analysis <- format(Sys.Date(), "%Y%b%d")
g.list <- vector(mode = "list", length = 5) # Create a plotting list container
var.i <- "BMI"
threshold.vec <- seq(1000, 5000, 1000)
for (threshold.i in threshold.vec) {
  i <- which(threshold.vec == threshold.i)
  cDEGs_num_threshold <- filter(cDEGs_num_all, threshold == threshold.i, variable == var.i)
  g <- cDEGs_num_threshold |>
    ggplot(aes(y = num.cDEGs, x = method, group = random.seed)) +
    geom_boxplot(outlier.shape = NA, mapping = aes(group = method), lwd = 0.8) +
    geom_line(aes(color = random.seed)) +
    geom_point(aes(color = random.seed), size = 0.8) + 
    # Annotate the max value of cDEGs num per seed by enlarging its point size
    geom_point(data = cDEGs_num_threshold |>
                 group_by(random.seed) |>
                 mutate(max.num = max(num.cDEGs)) |>
                 dplyr::filter(num.cDEGs == max.num),
               aes(y = max.num, x = method, color = random.seed), 
               size = 2) +
    facet_wrap(vars(threshold)) +
    labs(
      y = paste0("Number of cDEGs"),
      x = "Method",
      color = "Seed"
    ) +
    theme_BMA() +
    theme(
      axis.text.x = element_text(
        size = 9,
        angle = 90,
        hjust = 1,
        color = "black",
        vjust = 0.5
    ))
  
  if (threshold.i != 5000) {
    g.list[[i]] <- g
  } else {
    g_s <- g + theme(legend.position = c(0.2, 0.5))
    legend <- get_legend(g_s)
    g.list[[i]] <- g
    return(legend)
  }
}
g.final <- plot_grid(g.list[[1]], g.list[[2]], g.list[[3]], g.list[[4]], g.list[[5]], ncol = 3, labels = LETTERS[1:5], legend)
ggsave(
  filename = sprintf(
    "../ApplicationResult/AddViz/ratecDEGs/%s_%s_line_latex.eps",
    date.analysis, var.i
  ),
  plot = g.final,
  device = "eps",
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)


