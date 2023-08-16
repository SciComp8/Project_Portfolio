library(ggplot2)
library(ggpubr) 
library(viridis)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(readxl)
ensembl.id <- bitr("IRF3", 
                   fromType = "SYMBOL", 
                   toType = "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db")
ensembl.id
# SYMBOL         ENSEMBL
# 1   IRF3 ENSG00000126456

# Check the most frequent BMAseq best model across 10 seeds
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
# ~1+MHABNWBC

# IRF3 
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
dat.expr.train <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i))
dat.pheno.train <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i))
dat.expr.test <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i))
dat.pheno.test <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i))

# Use the whole population
dat.expr.all <- cbind(dat.expr.train, dat.expr.test)
dat.pheno.all <- rbind(dat.pheno.train, dat.pheno.test)
IRF3.expr.all <- dat.expr.all[rownames(dat.expr.all) == "ENSG00000126456.11", ]
IRF3.expr.all.t <- t(IRF3.expr.all)
IRF3.all <- data.frame(expr = IRF3.expr.all.t[, 1],
                       BMI = dat.pheno.all$BMI, 
                       AGE = dat.pheno.all$AGE,
                       SEX = dat.pheno.all$SEX,
                       MHABNWBC = dat.pheno.all$MHABNWBC)
IRF3.all$BMI_AGE_SEX_MHABNWBC <- paste(paste0("BMI", IRF3.all$BMI), 
                                       paste0("AGE", IRF3.all$AGE), 
                                       paste0("SEX", IRF3.all$SEX), 
                                       paste0("WBC", IRF3.all$MHABNWBC), sep = " + ")
IRF3.all$BMI_MHABNWBC <- paste(paste0("BMI", IRF3.all$BMI), 
                               paste0("WBC", IRF3.all$MHABNWBC), sep = " + ")
IRF3.all$AGE_MHABNWBC <- paste(paste0("AGE", IRF3.all$AGE), 
                               paste0("WBC", IRF3.all$MHABNWBC), sep = " + ")

# Set the plot theme
theme_BMA <- function(base_size = 12,
                      base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6), # Customize the box border in strip
      strip.placement = "inside"
    )
}

make_boxplot <- function(gene.symbol = "MTHFS", data.type = "train", var.x = "BMI", test.type = "wilcox.test") {
  p <- ggplot(data = get(paste0(gene.symbol, ".", data.type)),
              mapping = aes(x = get(var.x), y = expr)) +
    geom_boxplot(aes(group = get(var.x), color = get(var.x)),
                 outlier.shape = NA) + 
    geom_jitter(aes(color = get(var.x)), size = 0.8) + 
    stat_compare_means(method = test.type, 
                       label.x.npc = "middle",
                       # aes(label = ifelse(p < 1.e-3, "p < 0.001", sprintf("p = %4.3f", as.numeric(..p.format..))))
    ) + 
    theme_BMA() + 
    scale_color_viridis(discrete = T) + 
    ylab(paste0(gene.symbol, " expression level")) + 
    xlab(var.x)
  return(p)
}

p1 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "BMI", test.type = "wilcox.test")
p2 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "MHABNWBC", test.type = "wilcox.test")
p3 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "BMI_MHABNWBC", test.type = "kruskal.test")
p4 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "AGE_MHABNWBC", test.type = "kruskal.test")
p5 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "BMI_AGE_SEX_MHABNWBC", test.type = "kruskal.test")
g1 <- plot_grid(p1, p5, ncol = 2, labels = LETTERS[1:2])
g2 <- plot_grid(p2, p3, p4, ncol = 3, labels = LETTERS[3:5])

date.analysis <- format(Sys.Date(), "%Y%b%d")
ggsave(filename = sprintf("../ApplicationResult/AddViz/IRF3_BoxPlot/%s_%s_%s_%s.1.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g1, device = cairo_ps, dpi = 600, width = 10, height = 6, units = "in")

date.analysis <- format(Sys.Date(), "%Y%b%d")
ggsave(filename = sprintf("../ApplicationResult/AddViz/IRF3_BoxPlot/%s_%s_%s_%s.2.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g2, device = cairo_ps, dpi = 600, width = 12, height = 6, units = "in")
