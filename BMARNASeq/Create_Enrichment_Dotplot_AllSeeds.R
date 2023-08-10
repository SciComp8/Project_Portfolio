library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(latex2exp)

##------BMAseq-----
# Show all cDEGs uniquely identified by BMAseq that appear in at least one seed
var.name <- "BMI"
threshold.i <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
class_freq_all_seed_list <- vector(mode = "list", length = 10L)
names(class_freq_all_seed_list) <- seed.vec 
for (seed.i in seed.vec) {
  file.name <- paste0("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/", var.name, "_", threshold.i, "_", seed.i, ".RDS")
  class_freq_per_seed <- readRDS(file.name)
  class_freq_per_seed_BMAseq <- filter(class_freq_per_seed, Class == "100000000") |>
    dplyr::select(Class) |>
    mutate(Seed = seed.i) |>
    rownames_to_column("cDEG")
  class_freq_all_seed_list[[as.character(seed.i)]] <- class_freq_per_seed_BMAseq 
}

class_freq_all_seed_df <- do.call(rbind, class_freq_all_seed_list)

# Show all cDEGs uniquely identified by BMAseq that appear in at least 5 seeds
class_freq_all_seed_df_new <- summarize(.data = class_freq_all_seed_df, cDEG.freq = n(), .by = "cDEG")
class_freq_all_seed_df_select <- dplyr::filter(class_freq_all_seed_df_new, cDEG.freq >= 5)
unique_cDEG_all_seed <- sub("\\..*", "", class_freq_all_seed_df_select$cDEG)
unique_cDEG_all_seed_gene_symbol <- bitr(unique_cDEG_all_seed, fromType = "ENSEMBL", 
                                         toType = c("SYMBOL"), 
                                         OrgDb = "org.Hs.eg.db")

# Perform ORA analysis on cDEGs uniquely identified by BMAseq that appear in at least 5 seeds
ora.obj <-
  enrichGO(
    gene          = unique_cDEG_all_seed,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = T
  )

ora.obj.df <- ora.obj@result
ora.obj.ordered <- ora.obj.df[order(ora.obj.df$Count, decreasing = T), ]

p <- cnetplot(x = ora.obj,
         showCategory = ora.obj.ordered$Description[1:6],
         circular = T, 
         color.params = list(foldChange = NULL, edge = TRUE)) + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

date.analysis <- format(Sys.Date(), "%Y%b%d")
ggsave(filename = sprintf("../ApplicationResult/AddViz/cnetplot/%s_%s_%s_%s.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 11, height = 7, units = "in")

make_dotplot <- function(data.set = ora.obj.ordered, 
                         top.n = 20, 
                         method.name = "BMAseq", 
                         color.lower = NULL, 
                         color.upper = NULL, 
                         color.type = c("pvalue", "qscore")) {
  if (nrow(data.set) < top.n) {
    data.set <- data.set
  } else {
    data.set <- data.set[1:top.n, ]
  }
  term.ordered <- rev(data.set$Description)
  color.type <- match.arg(color.type, c("pvalue", "qscore"))
  min.count <- min(data.set$Count)
  max.count <- max(data.set$Count)
  
  if (is.null(color.lower)) {
    color.lower <- 0
  }
  if (is.null(color.upper)) {
    if (color.type == "pvalue") {
      color.upper <- range(data.set$pvalue)[2] |> round(2)
    } else if (color.type == "qscore") {
      data.set <- mutate(data.set, qscore = -log(pvalue, base = 10)) # Use the original p value to calculate the q score
      color.upper <- range(data.set$qscore)[2] |> round(2) + 0.5 # Increase the upper bound of q score
    }
  } 
  
  p <- ggplot(data = data.set, 
              mapping = aes(x = Count, y = Description)) + 
    geom_point(aes(size = Count, color = get(color.type))) +
    scale_colour_gradient(limits = c(color.lower, color.upper), low = "blue", high = "red") + 
    scale_size_continuous(breaks = seq(min.count, max.count, 1)) + # Count, not continuous scale
    scale_x_continuous(breaks = seq(min.count, max.count, 1)) + # Count, not continuous scale
    scale_y_discrete(limits = term.ordered, labels = scales::label_wrap(70)) + 
    labs(y = NULL,
         color = ifelse(color.type == "pvalue", "p-value", 
                        TeX("-$log_{10}$ (p-value)"))
         ) +
    ggtitle(paste0("Unique cDEGs of ", var.name)) + 
    theme_bw(base_size = 14, base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 13, color = "black")) + 
    guides(size = guide_legend(order = 1)) # Order the discrete legend
  
  return(p)
}


p <- make_dotplot(data.set = ora.obj.ordered, 
                  top.n = 20, 
                  method.name = "BMAseq", 
                  color.lower = NULL, 
                  color.upper = NULL, 
                  color.type = "qscore")

date.analysis <- format(Sys.Date(), "%Y%b%d")
ggsave(filename = sprintf("../ApplicationResult/AddViz/DotPlot/%s_%s_%s_%s.eps", date.analysis, "BMAseq", var.name, "all.seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 10, height = 6, units = "in")

##------edgeR_UVM-----
# Show all cDEGs uniquely identified by edgeR_UVM that appear in at least one seed
var.name <- "BMI"
threshold.i <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
class_freq_all_seed_list <- vector(mode = "list", length = 10L)
names(class_freq_all_seed_list) <- seed.vec 
for (seed.i in seed.vec) {
  ## Import the matrix that records the class of each cDEGs per seed
  file.name <- paste0("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/", var.name, "_", threshold.i, "_", seed.i, ".RDS")
  class_freq_per_seed <- readRDS(file.name)
  class_freq_per_seed_BMAseq <- filter(class_freq_per_seed, Class == "000100000") |>
    dplyr::select(Class) |>
    mutate(Seed = seed.i) |>
    rownames_to_column("cDEG")
  class_freq_all_seed_list[[as.character(seed.i)]] <- class_freq_per_seed_BMAseq 
}

class_freq_all_seed_df <- do.call(rbind, class_freq_all_seed_list)

# Show all cDEGs uniquely identified by edgeR_UVM that appear in at least 5 seeds 
class_freq_all_seed_df_new <- summarize(.data = class_freq_all_seed_df, cDEG.freq = n(), .by = "cDEG")
class_freq_all_seed_df_select <- dplyr::filter(class_freq_all_seed_df_new, cDEG.freq >= 5)
unique_cDEG_all_seed <- sub("\\..*", "", class_freq_all_seed_df_select$cDEG)
unique_cDEG_all_seed_gene_symbol <- bitr(unique_cDEG_all_seed, fromType = "ENSEMBL", 
                                         toType = c("SYMBOL"), 
                                         OrgDb = "org.Hs.eg.db")

# Perform ORA analysis on cDEGs uniquely identified by edgeR_UVM that appear in at least 5 seeds
ora.obj <-
  enrichGO(
    gene          = unique_cDEG_all_seed,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = T
  )

ora.obj.df <- ora.obj@result
ora.obj.ordered <- ora.obj.df[order(ora.obj.df$Count, decreasing = T), ]

p <- cnetplot(x = ora.obj,
         showCategory = ora.obj.ordered$Description[1:6],
         circular = T, 
         color.params = list(foldChange = NULL, edge = TRUE)) + 
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

date.analysis <- format(Sys.Date(), "%Y%b%d")
ggsave(filename = sprintf("../ApplicationResult/AddViz/cnetplot/%s_%s_%s_%s.eps", date.analysis, "edgeR_UVM", var.name, "all.seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 11, height = 7, units = "in")

p <- make_dotplot(data.set = ora.obj.ordered, 
                  top.n = 20, 
                  method.name = "edgeR_UVM", 
                  color.lower = NULL, 
                  color.upper = NULL, 
                  color.type = "qscore")

date.analysis <- format(Sys.Date(), "%Y%b%d")
ggsave(filename = sprintf("../ApplicationResult/AddViz/DotPlot/%s_%s_%s_%s.eps", date.analysis, "edgeR_UVM", var.name, "all.seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 10, height = 6, units = "in")

