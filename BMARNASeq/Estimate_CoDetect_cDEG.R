var.name <- c("BMI")

easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", 
                        "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", 
                        "openxlsx", "grDevices", "viridis", "cowplot", "clusterProfiler",
                        "org.Hs.eg.db", "latex2exp", "enrichplot") |> suppressPackageStartupMessages()

# Estimate the number of genes co-identified by one approach and at least one of the other approaches - one seed
seed <- 8809678
class.freq.8809678 <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_%s.RDS", var.name, seed))
class.freq.8809678 <- class.freq.8809678[, -ncol(class.freq.8809678)]
class.freq.8809678 <- class.freq.8809678 |> mutate(across(.cols = everything(), as.numeric))
# class.freq.8809678[nrow(class.freq.8809678) + 1, ] <- colSums(class.freq.8809678)
# rownames(class.freq.8809678)[nrow(class.freq.8809678)] <- "col.sum"
class.freq.8809678[, ncol(class.freq.8809678) + 1] <- rowSums(class.freq.8809678)
colnames(class.freq.8809678)[ncol(class.freq.8809678)] <- "row.sum"


r.1 <- class.freq.8809678 |> 
  filter(row.sum >= 2) |> 
  dplyr::select(!row.sum) |>
  colSums()
r.1 <- c(r.1, seed = 8809678)

r.2 <- class.freq.8809678 |> 
  filter(row.sum >= 2) |> 
  dplyr::select(!row.sum) |>
  colSums()
r.2 <- c(r.2, seed = 8809678)

rbind(r.1, r.2)

# Estimate the number of genes co-identified by one approach and at least one of the other approaches - all seeds
r.mat <- matrix(0, nrow = 10, ncol = 10)
colnames(r.mat) <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM",
                     "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM", "Seed")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)

for (i in 1:10) {
  class.freq <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_%s.RDS", var.name, seed.vec[i]))
  class.freq <- class.freq[, -ncol(class.freq)]
  class.freq <- class.freq |> mutate(across(.cols = everything(), as.numeric))
  class.freq[, ncol(class.freq) + 1] <- rowSums(class.freq)
  colnames(class.freq)[ncol(class.freq)] <- "row.sum"
  
  r <- class.freq |> 
    filter(row.sum >= 2) |> 
    dplyr::select(!row.sum) |>
    colSums()
  r <- c(r, seed = seed.vec[i])
  
  r.mat[i, ] <- r
}

r.df <- data.frame(r.mat) |>
  pivot_longer(cols = !Seed, names_to = "Method", values_to = "Number")

r.df$Method <- factor(r.df$Method,
                      levels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM"))

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

date.analysis <- format(Sys.Date(), "%Y%b%d")
p <- r.df |>
  ggplot(mapping = aes(y = Number, x = Method, color = Method)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  geom_jitter(aes(color = Method), size = 0.8) + 
  labs(y = "Number of cDEGs co-detected with at least 2 approaches", x = "Method") + 
  theme_BMA() + 
  scale_color_viridis(discrete = T)
ggsave(filename = sprintf("../ApplicationResult/AddViz/shared_discovery/%s_%s_%s.eps", date.analysis, var.name, "all_seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 7, height = 8, units = "in")


# Estimate the number of genes co-identified by one approach and a certain number of the other approaches - all seeds
r_fun <- function(num.method = NULL) {
  r <- class.freq |> 
    filter(row.sum == num.method) |> 
    dplyr::select(!row.sum) |>
    colSums()
  r <- c(r, seed = seed.vec[i], num.method = num.method)
}

r.list <- vector(mode = "list", length = 8)
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)

for (z in 1:length(r.list)) {
  r.mat <- matrix(0, nrow = 10, ncol = 11)
  colnames(r.mat) <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM",
                       "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM", "Seed", "Number_method")
  for (i in 1:10) {
    class.freq <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_%s.RDS", var.name, seed.vec[i]))
    class.freq <- class.freq[, -ncol(class.freq)]
    class.freq <- class.freq |> mutate(across(.cols = everything(), as.numeric))
    class.freq[, ncol(class.freq) + 1] <- rowSums(class.freq)
    colnames(class.freq)[ncol(class.freq)] <- "row.sum"
    
    r <- r_fun(num.method = z + 1)
    
    r.mat[i, ] <- r
  }
  r.list[[z]] <- r.mat
}

r.list.to.df <- do.call(rbind, r.list) |> 
  data.frame()
r.list.to.df.long <- r.list.to.df |> 
  pivot_longer(cols = !c(Seed, Number_method), names_to = "Method", values_to = "Number")

r.list.to.df.long$Method <- factor(r.list.to.df.long$Method,
                                   levels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM"))

date.analysis <- format(Sys.Date(), "%Y%b%d")
p <- r.list.to.df.long |>
  ggplot(mapping = aes(y = Number, x = Method, color = Method)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  geom_jitter(aes(color = Method), size = 0.8) + 
  facet_wrap(~Number_method) + 
  labs(y = "Number of cDEGs co-detected with varying numbers of approaches", x = "Method") + 
  theme_BMA() + 
  scale_color_viridis(discrete = T)
ggsave(filename = sprintf("../ApplicationResult/AddViz/shared_discovery/%s_%s_%s.eps", date.analysis, var.name, "vary_no_method_all_seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 10, height = 8, units = "in")


# Estimate the number of genes uniquely identified by one approach - all seeds
r.mat <- matrix(0, nrow = 10, ncol = 10)
colnames(r.mat) <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM",
                     "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM", "Seed")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)

for (i in 1:10) {
  class.freq <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_%s.RDS", var.name, seed.vec[i]))
  class.freq <- class.freq[, -ncol(class.freq)]
  class.freq <- class.freq |> mutate(across(.cols = everything(), as.numeric))
  class.freq[, ncol(class.freq) + 1] <- rowSums(class.freq)
  colnames(class.freq)[ncol(class.freq)] <- "row.sum"
  
  r <- class.freq |> 
    filter(row.sum == 1) |> 
    dplyr::select(!row.sum) |>
    colSums()
  r <- c(r, seed = seed.vec[i])
  
  r.mat[i, ] <- r
}

r.df <- data.frame(r.mat) |>
  pivot_longer(cols = !Seed, names_to = "Method", values_to = "Number")

r.df$Method <- factor(r.df$Method,
                      levels = c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM"))

date.analysis <- format(Sys.Date(), "%Y%b%d")
p <- r.df |>
  ggplot(mapping = aes(y = Number, x = Method, color = Method)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  geom_jitter(aes(color = Method), size = 0.8) + 
  labs(y = "Number of cDEGs uniquely detected with an approach", x = "Method") + 
  theme_BMA() + 
  scale_color_viridis(discrete = T)
ggsave(filename = sprintf("../ApplicationResult/AddViz/unique_discovery/%s_%s_%s.eps", date.analysis, var.name, "all_seed"),
       plot = p, device = cairo_ps, dpi = 600, width = 6, height = 7, units = "in")
