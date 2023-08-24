##------Attach libraries------
library(ggplot2)
library(ggpubr) 
library(viridis)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(readxl)
options(scipen = 999) # Remove scientific notation

##------Map gene symbols to ENSEMBL IDs------
ensembl.id <- bitr("CXCL10", 
                   fromType = "SYMBOL", 
                   toType = "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db")
ensembl.id
# SYMBOL         ENSEMBL
# 1   CXCL10 ENSG00000169245

##------Check the most frequent BMAseq best model across 10 seeds------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")

##------Use the whole data to make the boxplot------
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.i <- 120
dat.expr.train <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i))
dat.pheno.train <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i))
dat.expr.test <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i))
dat.pheno.test <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i))
dat.expr.all <- cbind(dat.expr.train, dat.expr.test)
dat.pheno.all <- rbind(dat.pheno.train, dat.pheno.test)
CXCL10.expr.all <- dat.expr.all[grepl(pattern = "ENSG00000169245", x = rownames(dat.expr.all)), ]
CXCL10.expr.all.t <- t(CXCL10.expr.all)
CXCL10.all <- data.frame(expr = CXCL10.expr.all.t[, 1],
                         BMI = dat.pheno.all$BMI, 
                         AGE = dat.pheno.all$AGE,
                         SEX = dat.pheno.all$SEX,
                         MHABNWBC = dat.pheno.all$MHABNWBC)

##------Use the data from DESeq2 dataset to make the boxplot------
cts.all <- as.matrix(dat.expr.all)
coldata.all <- dat.pheno.all

# Transform the TMM normalization factors to be used in DESeq2
lib.size <- colSums(cts.all)
norm.factor <- calcNormFactors(cts.all, method = "TMM")
size.factor <- lib.size*norm.factor/exp(mean(log(lib.size*norm.factor))) 

var <- "BMI"
dds <- DESeqDataSetFromMatrix(countData = cts.all, 
                              colData = coldata.all, 
                              design = formula(paste0("~", var)))
sizeFactors(dds) <- size.factor
results(DESeq(dds))
dds_vst <- vst(dds)

dds_vst_df <- assay(dds_vst) |> data.frame()
CXCL10.expr.all.deseq2 <- assay(dds_vst)[grepl(pattern = "ENSG00000169245", x = rownames(assay(dds_vst))), ]
CXCL10.all.deseq2 <- data.frame(expr = CXCL10.expr.all.deseq2,
                                BMI = dat.pheno.all$BMI, 
                                AGE = dat.pheno.all$AGE,
                                SEX = dat.pheno.all$SEX,
                                MHABNWBC = dat.pheno.all$MHABNWBC)

##------Set the factor levels of phenotypes------
CXCL10.all$BMI <- factor(CXCL10.all$BMI, levels = c("low", "high"))
CXCL10.all$AGE <- factor(CXCL10.all$AGE, levels = c("young", "old"))
CXCL10.all$SEX <- factor(CXCL10.all$SEX, levels = c("female", "male"))
CXCL10.all$MHABNWBC <- factor(CXCL10.all$MHABNWBC, levels = c("no", "yes"))

CXCL10.all$BMI_AGE_SEX_MHABNWBC <- paste(paste0("AGE", CXCL10.all$AGE), 
                                         paste0("SEX", CXCL10.all$SEX), 
                                         paste0("WBC", CXCL10.all$MHABNWBC),
                                         paste0("BMI", CXCL10.all$BMI), sep = " + ")
combs <- table(CXCL10.all$BMI_AGE_SEX_MHABNWBC) |> names()
combs <- combs[c(9:16, 1:8)] 
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
CXCL10.all$BMI_AGE_SEX_MHABNWBC <- factor(CXCL10.all$BMI_AGE_SEX_MHABNWBC,
                                          levels = combs_new)


CXCL10.all$BMI_MHABNWBC <- paste(paste0("WBC", CXCL10.all$MHABNWBC), 
                                 paste0("BMI", CXCL10.all$BMI), sep = " + ")
combs <- table(CXCL10.all$BMI_MHABNWBC) |> names()
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
CXCL10.all$BMI_MHABNWBC <- factor(CXCL10.all$BMI_MHABNWBC,
                                  levels = combs_new)


CXCL10.all$AGE_MHABNWBC <- paste(paste0("AGE", CXCL10.all$AGE), 
                                 paste0("WBC", CXCL10.all$MHABNWBC), sep = " + ")
combs <- table(CXCL10.all$AGE_MHABNWBC) |> names()
combs <- combs[c(3:4, 1:2)] 
CXCL10.all$AGE_MHABNWBC <- factor(CXCL10.all$AGE_MHABNWBC,
                                levels = combs)

CXCL10.all$BMI_AGE_MHABNWBC <- paste(paste0("AGE", CXCL10.all$AGE), 
                                     paste0("WBC", CXCL10.all$MHABNWBC),
                                     paste0("BMI", CXCL10.all$BMI), sep = " + ")
combs <- table(CXCL10.all$BMI_AGE_MHABNWBC) |> names()
combs <- combs[c(5:8, 1:4)] 
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
CXCL10.all$BMI_AGE_MHABNWBC <- factor(CXCL10.all$BMI_AGE_MHABNWBC,
                                      levels = combs_new)


##------Set the factor levels of phenotypes in vst transformed dataset------
CXCL10.all.deseq2$BMI <- factor(CXCL10.all.deseq2$BMI, levels = c("low", "high"))
CXCL10.all.deseq2$AGE <- factor(CXCL10.all.deseq2$AGE, levels = c("young", "old"))
CXCL10.all.deseq2$SEX <- factor(CXCL10.all.deseq2$SEX, levels = c("female", "male"))
CXCL10.all.deseq2$MHABNWBC <- factor(CXCL10.all.deseq2$MHABNWBC, levels = c("no", "yes"))

CXCL10.all.deseq2$BMI_AGE_SEX_MHABNWBC <- paste(paste0("AGE", CXCL10.all.deseq2$AGE), 
                                                paste0("SEX", CXCL10.all.deseq2$SEX), 
                                                paste0("WBC", CXCL10.all.deseq2$MHABNWBC),
                                                paste0("BMI", CXCL10.all.deseq2$BMI), sep = " + ")
combs <- table(CXCL10.all.deseq2$BMI_AGE_SEX_MHABNWBC) |> names()
combs <- combs[c(9:16, 1:8)] 
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
CXCL10.all.deseq2$BMI_AGE_SEX_MHABNWBC <- factor(CXCL10.all.deseq2$BMI_AGE_SEX_MHABNWBC,
                                               levels = combs_new)


CXCL10.all.deseq2$BMI_MHABNWBC <- paste(paste0("WBC", CXCL10.all.deseq2$MHABNWBC), 
                                        paste0("BMI", CXCL10.all.deseq2$BMI), sep = " + ")
combs <- table(CXCL10.all.deseq2$BMI_MHABNWBC) |> names()
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
CXCL10.all.deseq2$BMI_MHABNWBC <- factor(CXCL10.all.deseq2$BMI_MHABNWBC,
                                         levels = combs_new)


CXCL10.all.deseq2$AGE_MHABNWBC <- paste(paste0("AGE", CXCL10.all.deseq2$AGE), 
                                      paste0("WBC", CXCL10.all.deseq2$MHABNWBC), sep = " + ")
combs <- table(CXCL10.all.deseq2$AGE_MHABNWBC) |> names()
combs <- combs[c(3:4, 1:2)] 
CXCL10.all.deseq2$AGE_MHABNWBC <- factor(CXCL10.all.deseq2$AGE_MHABNWBC,
                                       levels = combs)

CXCL10.all.deseq2$BMI_AGE_MHABNWBC <- paste(paste0("AGE", CXCL10.all.deseq2$AGE), 
                                            paste0("WBC", CXCL10.all.deseq2$MHABNWBC),
                                            paste0("BMI", CXCL10.all.deseq2$BMI), sep = " + ")
combs <- table(CXCL10.all.deseq2$BMI_AGE_MHABNWBC) |> names()
combs <- combs[c(5:8, 1:4)] 
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
CXCL10.all.deseq2$BMI_AGE_MHABNWBC <- factor(CXCL10.all.deseq2$BMI_AGE_MHABNWBC,
                                             levels = combs_new)

##------Customize the plot theme------
theme_BMA <- function(base_size = 14,
                      base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, color = "black", vjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6), # Customize the box border in strip
      strip.placement = "inside"
    )
}

##------Make a boxplot function------
make_boxplot <- function(gene.symbol = "MTHFS", 
                         data.type = "train", 
                         var.x = "BMI", 
                         test.type = "wilcox.test",
                         paired_comparisons_list = NULL,
                         y_transform = c("origin", "log2", "vst")) {
  
  data_set <- get(paste0(gene.symbol, ".", data.type))
  data_set[["expr"]] <- data_set[["expr"]] + 0.5
  
  y_transform <- match.arg(y_transform, c("origin", "log2", "vst"))
  if (y_transform == "origin") {
    y_lab_name <- paste0(gene.symbol, " count + 0.5")
  } else if (y_transform == "log2") {
    data_set[["expr"]] <- log2(data_set[["expr"]]) # Log2 transformed gene expression
    y_lab_name <- bquote(paste("log"[2]*" (", .(gene.symbol), " count + 0.5)"))
  } else if (y_transform == "vst") {
    data_set <- get(paste0(gene.symbol, ".", data.type, ".deseq2"))
    y_lab_name <- paste0(gene.symbol, " normalized count")
  }
  
  p <- ggplot(data = data_set,
              mapping = aes(x = get(var.x), y = expr)) +
    geom_boxplot(aes(group = get(var.x), color = get(var.x)),
                 lwd = 1.2,
                 outlier.shape = NA) + 
    geom_jitter(aes(color = get(var.x)), size = 0.8) + 
    stat_compare_means(method = test.type, # This is the overall comparison
                       label.x.npc = "middle",
                       size = 5) + 
    theme_BMA() + 
    scale_color_viridis(discrete = T) + 
    xlab(gsub(pattern = "_", replacement = " + ", x = var.x)) + 
    ylab(y_lab_name)
  
  if (is.null(paired_comparisons_list)) {
    p
  } else {
    p <- p + stat_compare_means(label.x.npc = "middle",
                                comparisons = paired_comparisons_list,
                                method = "wilcox.test")
  }
}

##------Make boxplots------
p1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI", test.type = "wilcox.test", y_transform = "log2")
p1.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI", test.type = "wilcox.test", y_transform = "vst")
p2 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "MHABNWBC", test.type = "wilcox.test", y_transform = "log2")
p2.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "MHABNWBC", test.type = "wilcox.test", y_transform = "vst")
p3 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI_MHABNWBC", 
                   test.type = "kruskal.test", y_transform = "log2")
p3.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI_MHABNWBC", 
                     test.type = "kruskal.test", y_transform = "vst")
p4 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "AGE_MHABNWBC", test.type = "kruskal.test", y_transform = "log2")
p4.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "AGE_MHABNWBC", test.type = "kruskal.test", y_transform = "vst")
p5 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI_AGE_MHABNWBC", 
                   test.type = "kruskal.test", y_transform = "log2")
p5.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI_AGE_MHABNWBC", 
                     test.type = "kruskal.test", y_transform = "vst")
p6 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI_AGE_SEX_MHABNWBC", test.type = "kruskal.test", y_transform = "log2")
p6.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "BMI_AGE_SEX_MHABNWBC", test.type = "kruskal.test", y_transform = "vst")
p7 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "AGE", test.type = "wilcox.test", y_transform = "log2")
p7.1 <- make_boxplot(gene.symbol = "CXCL10", data.type = "all", var.x = "AGE", test.type = "wilcox.test", y_transform = "vst")

##------Put plots together in a panel------
g1 <- plot_grid(p1, p6, ncol = 2)
g1.1 <- plot_grid(p1.1, p6.1, ncol = 2)
g2 <- plot_grid(p2, p5, ncol = 2)
g2.1 <- plot_grid(p2.1, p5.1, ncol = 2)
g3 <- plot_grid(p7, p4, ncol = 2)
g3.1 <- plot_grid(p7.1, p4.1, ncol = 2)

##------Save plots------
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.name <- "BMI"
ggsave(filename = sprintf("../ApplicationResult/AddViz/CXCL10_BoxPlot/%s_%s_%s_%s.1.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g1, device = cairo_ps, dpi = 600, width = 14, height = 7, units = "in")
ggsave(filename = sprintf("../ApplicationResult/AddViz/CXCL10_BoxPlot/%s_%s_%s_%s.1.1.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g1.1, device = cairo_ps, dpi = 600, width = 14, height = 7, units = "in")

ggsave(filename = sprintf("../ApplicationResult/AddViz/CXCL10_BoxPlot/%s_%s_%s_%s.2.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g2, device = cairo_ps, dpi = 600, width = 14, height = 7, units = "in")
ggsave(filename = sprintf("../ApplicationResult/AddViz/CXCL10_BoxPlot/%s_%s_%s_%s.2.1.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g2.1, device = cairo_ps, dpi = 600, width = 14, height = 7, units = "in")

ggsave(filename = sprintf("../ApplicationResult/AddViz/CXCL10_BoxPlot/%s_%s_%s_%s.3.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g3, device = cairo_ps, dpi = 600, width = 12, height = 6, units = "in")
ggsave(filename = sprintf("../ApplicationResult/AddViz/CXCL10_BoxPlot/%s_%s_%s_%s.3.1.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g3.1, device = cairo_ps, dpi = 600, width = 14, height = 6, units = "in")
