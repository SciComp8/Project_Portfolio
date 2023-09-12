image.date <- format(Sys.Date(), "%Y%b%d")
save.image(file = paste0(image.date, "IRF3_boxplot_ALiu.RData"))

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
ensembl.id <- bitr("IRF3", 
                   fromType = "SYMBOL", 
                   toType = "ENSEMBL", 
                   OrgDb = "org.Hs.eg.db")
ensembl.id
# SYMBOL         ENSEMBL
# 1   IRF3 ENSG00000126456

##------Check the most frequent BMAseq best model across 10 seeds------
best.model <- read_excel("../ApplicationData/derived/RandomSeed/BestModelMatrix/BMI_5000.xlsx")
# ~1+MHABNWBC

##------Use the whole count data to make the boxplot------
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.i <- 120
dat.expr.train <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i))
dat.pheno.train <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i))
dat.expr.test <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i))
dat.pheno.test <- readRDS(sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i))
dat.expr.all <- cbind(dat.expr.train, dat.expr.test)
dat.pheno.all <- rbind(dat.pheno.train, dat.pheno.test)
IRF3.expr.all <- dat.expr.all[rownames(dat.expr.all) == "ENSG00000126456.11", ]
IRF3.expr.all.t <- t(IRF3.expr.all)
IRF3.all <- data.frame(expr = IRF3.expr.all.t[, 1],
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

# assay(dds_vst) |> View() # View the transformed count
# b@assays@data@listData[[1]] 

dds_vst_df <- assay(dds_vst) |> data.frame()
IRF3.expr.all.deseq2 <- assay(dds_vst)[rownames(assay(dds_vst)) == "ENSG00000126456.11", ]
IRF3.all.deseq2 <- data.frame(expr = IRF3.expr.all.deseq2,
                              BMI = dat.pheno.all$BMI, 
                              AGE = dat.pheno.all$AGE,
                              SEX = dat.pheno.all$SEX,
                              MHABNWBC = dat.pheno.all$MHABNWBC)

##------Set the factor levels of phenotypes in count dataset------
IRF3.all$BMI <- factor(IRF3.all$BMI, levels = c("low", "high"))
IRF3.all$AGE <- factor(IRF3.all$AGE, levels = c("young", "old"))
IRF3.all$SEX <- factor(IRF3.all$SEX, levels = c("female", "male"))
IRF3.all$MHABNWBC <- factor(IRF3.all$MHABNWBC, levels = c("no", "yes"))


IRF3.all$BMI_AGE_SEX_MHABNWBC <- paste(paste0("AGE", IRF3.all$AGE), 
                                       paste0("SEX", IRF3.all$SEX), 
                                       paste0("WBC", IRF3.all$MHABNWBC),
                                       paste0("BMI", IRF3.all$BMI), sep = " + ")
combs <- table(IRF3.all$BMI_AGE_SEX_MHABNWBC) |> names()
combs <- combs[c(9:16, 1:8)] 
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
IRF3.all$BMI_AGE_SEX_MHABNWBC <- factor(IRF3.all$BMI_AGE_SEX_MHABNWBC,
                                        levels = combs_new)


IRF3.all$BMI_MHABNWBC <- paste(paste0("WBC", IRF3.all$MHABNWBC), 
                               paste0("BMI", IRF3.all$BMI), sep = " + ")
combs <- table(IRF3.all$BMI_MHABNWBC) |> names()
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
IRF3.all$BMI_MHABNWBC <- factor(IRF3.all$BMI_MHABNWBC,
                                levels = combs_new)


IRF3.all$AGE_MHABNWBC <- paste(paste0("AGE", IRF3.all$AGE), 
                               paste0("WBC", IRF3.all$MHABNWBC), sep = " + ")
combs <- table(IRF3.all$AGE_MHABNWBC) |> names()
IRF3.all$AGE_MHABNWBC <- factor(IRF3.all$AGE_MHABNWBC,
                                levels = c(combs[grep(pattern = "WBCno", x = combs)],
                                           combs[grep(pattern = "WBCyes", x = combs)]))


##------Set the factor levels of phenotypes in vst transformed dataset------
IRF3.all.deseq2$BMI <- factor(IRF3.all.deseq2$BMI, levels = c("low", "high"))
IRF3.all.deseq2$AGE <- factor(IRF3.all.deseq2$AGE, levels = c("young", "old"))
IRF3.all.deseq2$SEX <- factor(IRF3.all.deseq2$SEX, levels = c("female", "male"))
IRF3.all.deseq2$MHABNWBC <- factor(IRF3.all.deseq2$MHABNWBC, levels = c("no", "yes"))

IRF3.all.deseq2$BMI_AGE_SEX_MHABNWBC <- paste(paste0("AGE", IRF3.all.deseq2$AGE), 
                                              paste0("SEX", IRF3.all.deseq2$SEX), 
                                              paste0("WBC", IRF3.all.deseq2$MHABNWBC),
                                              paste0("BMI", IRF3.all.deseq2$BMI), sep = " + ")
combs <- table(IRF3.all.deseq2$BMI_AGE_SEX_MHABNWBC) |> names()
combs <- combs[c(9:16, 1:8)] 
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
IRF3.all.deseq2$BMI_AGE_SEX_MHABNWBC <- factor(IRF3.all.deseq2$BMI_AGE_SEX_MHABNWBC,
                                               levels = combs_new)

IRF3.all.deseq2$BMI_MHABNWBC <- paste(paste0("WBC", IRF3.all.deseq2$MHABNWBC), 
                                      paste0("BMI", IRF3.all.deseq2$BMI), sep = " + ")
combs <- table(IRF3.all.deseq2$BMI_MHABNWBC) |> names()
combs_new <- c()
combs_new[seq(1, length(combs), by = 2)] <- combs[c(F, T)]
combs_new[seq(2, length(combs), by = 2)] <- combs[c(T, F)]
IRF3.all.deseq2$BMI_MHABNWBC <- factor(IRF3.all.deseq2$BMI_MHABNWBC,
                                       levels = combs_new)


IRF3.all.deseq2$AGE_MHABNWBC <- paste(paste0("AGE", IRF3.all.deseq2$AGE), 
                                      paste0("WBC", IRF3.all.deseq2$MHABNWBC), sep = " + ")
combs <- table(IRF3.all.deseq2$AGE_MHABNWBC) |> names()
IRF3.all.deseq2$AGE_MHABNWBC <- factor(IRF3.all.deseq2$AGE_MHABNWBC,
                                       levels = c(combs[grep(pattern = "WBCno", x = combs)],
                                                  combs[grep(pattern = "WBCyes", x = combs)]))


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
    y_lab_name <- paste0("*", gene.symbol, "* ", " transformed count") # Italicize the gene symbol
  }
  
  p <- ggplot(data = data_set,
              mapping = aes(x = get(var.x), y = expr)) +
    geom_boxplot(aes(group = get(var.x), color = get(var.x)),
                 lwd = 1,
                 outlier.shape = NA) + 
    geom_jitter(aes(color = get(var.x)), size = 0.8) + 
    stat_compare_means(method = test.type, # This is the overall comparison
                       label.x.npc = 0.2,
                       size = 4) + 
    theme_BMA() + 
    theme(axis.title.y = ggtext::element_markdown()) + 
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
p1 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "BMI", test.type = "wilcox.test", y_transform = "vst")
p2 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "BMI_AGE_SEX_MHABNWBC", test.type = "kruskal.test", y_transform = "vst")
p3 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "MHABNWBC", test.type = "wilcox.test", y_transform = "vst")
p4 <- make_boxplot(gene.symbol = "IRF3", data.type = "all", var.x = "BMI_MHABNWBC", 
                   test.type = "kruskal.test", y_transform = "vst")
g1 <- plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])
g2 <- plot_grid(p3, p4, ncol = 2, labels = LETTERS[3:4])

##------Save plots------
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.name <- "BMI"
ggsave(filename = sprintf("../ApplicationResult/AddViz/IRF3_BoxPlot/%s_%s_%s_%s.latex.1.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g1, device = cairo_ps, dpi = 600, width = 8, height = 6, units = "in")
ggsave(filename = sprintf("../ApplicationResult/AddViz/IRF3_BoxPlot/%s_%s_%s_%s.latex.2.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = g2, device = cairo_ps, dpi = 600, width = 8, height = 6, units = "in")
