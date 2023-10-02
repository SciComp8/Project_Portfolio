##------Set up the environment------
box::use(dplyr = dplyr)
library(ggplot2)
source("theme_BMA.R")
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.name <- "BMI"

##------Preprocess data for mean rate of BMI-related cDEGs------
cDEGs_rate_uni <- readRDS("../ApplicationResult/Uni/RandomSeed/ratecDEGs/Boxplot/boxplot.data.rate.BMA.RDS")
cDEGs_rate_uni <- cDEGs_rate_uni[cDEGs_rate_uni$method != "BMAseq", ]
cDEGs_rate_multi <- readRDS("../ApplicationResult/Multi/RandomSeed/ratecDEGs/Boxplot/boxplot.data.rate.RDS")
cDEGs_rate_uni$method <- paste0(cDEGs_rate_uni$method, "_UVM")
cDEGs_rate_multi$method <- ifelse(cDEGs_rate_multi$method != "BMAseq", paste0(cDEGs_rate_multi$method, "_MVM"), cDEGs_rate_multi$method)

cDEGs_rate_all <- rbind(cDEGs_rate_uni, cDEGs_rate_multi)
cDEGs_rate_all$method <- factor(cDEGs_rate_all$method,
                                levels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM"))
cDEGs_rate_all$variable <- factor(cDEGs_rate_all$variable,levels = c("BMI","AGE","SEX","MHABNWBC"))
cDEGs_rate_all$random.seed <- factor(cDEGs_rate_all$random.seed,
                                     levels = c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999))
cDEGs_rate_BMI <- dplyr$filter(.data = cDEGs_rate_all, variable == "BMI")
cDEGs_rate_mean_BMI <- dplyr$summarise(.data = cDEGs_rate_BMI,
                                       mean.rate.cDEGs = mean(rate.cDEGs),
                                       .by = c(threshold, method))
cDEGs_rate_mean_BMI$threshold <- factor(cDEGs_rate_mean_BMI$threshold,
                                        levels = 1000:5000)
cDEGs_rate_mean_BMI$line.type <- ifelse(grepl(pattern = "BMAseq|*.UVM", x = cDEGs_rate_mean_BMI$method), "solid", "dashed") |> factor(levels = c("solid", "dashed"))

##------Make cDEGs mean rate line plot------
g1 <- cDEGs_rate_mean_BMI |>
  ggplot(aes(y = mean.rate.cDEGs, x = threshold, color = method, group = method, linetype = line.type)) +
  geom_line() +
  geom_point(size = 0.8) +
  labs(y = paste0("Mean rate of cDEGs"),
       x = "Ranking threshold",
       color = "Method",
       linetype = "Method") +
  theme_BMA() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      size = 9,
      angle = 0,
      hjust = 0.5,
      color = "black",
      vjust = 0.5
    )) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     limits = c(0.1, 0.4), # Modify the range
                     breaks = seq(0.1, 0.4, 0.05)) +
  scale_color_manual(values = c("#A73030FF", rep(c("#003C67FF", "#0073C2FF", "#8F7700FF", "#EFC000FF"), each = 2))) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", rep(c("solid", "dashed"), times = 4))))) + 
  scale_linetype(guide = "none")

##------Preprocess data for mean rediscovery rate of BMI-related cDEGs------
duplicate_rate_uni <- readRDS("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/UniModel/organized_duplicate_rate_data.noBMA.RDS")
duplicate_rate_multi <- readRDS("../ApplicationData/derived/RandomSeed/DuplicatedRateMatrix/MultiModel/organized_duplicate_rate_data.RDS")

duplicate_rate_uni$method <- paste0(duplicate_rate_uni$method, "_UVM")
duplicate_rate_multi$method <- ifelse(duplicate_rate_multi$method != "BMAseq", paste0(duplicate_rate_multi$method, "_MVM"), duplicate_rate_multi$method)
duplicate_rate_all <- rbind(duplicate_rate_uni, duplicate_rate_multi)
duplicate_rate_all$method <- factor(duplicate_rate_all$method,
                                    levels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "VoomLimma_UVM", "VoomLimma_MVM"),
                                    labels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM"))
duplicate_rate_all$variable <- factor(duplicate_rate_all$variable,levels = c("BMI","AGE","SEX","MHABNWBC"))

##------Make cDEGs mean rediscovery rate line and box plot------
duplicate_rate_all$threshold <- factor(duplicate_rate_all$threshold, levels = c(1000, 2000, 3000, 4000, 5000))
duplicate_rate_bmi <- dplyr$filter(duplicate_rate_all, variable == "BMI")
g2 <- duplicate_rate_bmi |>
  ggplot(mapping = aes(y = cDEGs.duplicate.rate.mean, x = method, group = threshold)) + 
  geom_boxplot(outlier.shape = NA, mapping = aes(group = method), lwd = 0.8) +
  geom_line(aes(color = threshold)) + 
  geom_point(aes(color = threshold), size = 0.8) + 
  # Annotate the max value of cDEGs num per seed by enlarging its point size
  geom_point(
    data = duplicate_rate_bmi |>
      dplyr$group_by(threshold) |>
      dplyr$mutate(max.rate.mean = max(cDEGs.duplicate.rate.mean)) |>
      dplyr$filter(cDEGs.duplicate.rate.mean == max.rate.mean),
    aes(y = max.rate.mean, x = method, color = threshold),
    size = 2
  ) +
  coord_cartesian(clip = "off") + 
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     breaks = seq(0, 1, 0.05)) + 
  theme_BMA() + 
  scale_color_manual(values = c("#A20056FF", "#008B45FF", "#008280FF", "#3B4992FF", "#631879FF")) + # Avoid red-green combinations difficult for color-blind readers 
  labs(y = "Mean rediscovery rate of cDEGs", x = "Method", color = "Ranking threshold") + 
  theme(legend.position = "right")

ggsave(filename = sprintf("../ApplicationResult/AddViz/2_MeanRates/mean_rate_%s_%s_%s.latex.eps", date.analysis, "BMAseq", var.name),
       plot = g1, device = cairo_ps, dpi = 600, width = 5, height = 4, units = "in")
ggsave(filename = sprintf("../ApplicationResult/AddViz/2_MeanRates/mean_rediscovery_rate_%s_%s_%s.latex.eps", date.analysis, "BMAseq", var.name),
       plot = g2, device = cairo_ps, dpi = 600, width = 5, height = 4, units = "in")
