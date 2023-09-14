##------Set up the environment------
var.name <- c("BMI")
easypackages::libraries("tidyverse", "grDevices", "viridis", "cowplot") |> suppressPackageStartupMessages()
source("theme_BMA.R")
date.analysis <- format(Sys.Date(), "%Y%b%d")

##------Estimate the number of genes co-identified by a particular approach and at least one of the other approaches across all random seed trials------
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

p.shared.discovery <- r.df |>
  ggplot(mapping = aes(y = Number, x = Method, color = Method)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  geom_jitter(aes(color = Method), size = 0.8) + 
  labs(y = "Number of cDEGs co-detected with at least 2 methods", x = "Method") + 
  theme_BMA() + 
  scale_color_viridis(discrete = T)

##------Estimate the number of genes uniquely identified by one approach across all random seed trials------
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

p.unique.discovery <- r.df |>
  ggplot(mapping = aes(y = Number, x = Method, color = Method)) + 
  geom_boxplot(outlier.shape = NA, lwd = 1.2) + 
  geom_jitter(aes(color = Method), size = 0.8) + 
  labs(y = "Number of cDEGs uniquely detected with a method", x = "Method") + 
  theme_BMA() + 
  scale_color_viridis(discrete = T)

##------Output the final figure panel------
g <- plot_grid(p.shared.discovery, p.unique.discovery, ncol = 2, labels = LETTERS[1:2])
ggsave(filename = sprintf("../ApplicationResult/AddViz/unique_discovery/%s_%s_%s.eps", date.analysis, var.name, "all_seed"),
       plot = g, device = cairo_ps, dpi = 600, width = 8, height = 6, units = "in")
