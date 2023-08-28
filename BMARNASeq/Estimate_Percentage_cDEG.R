#####Among cDEGs detected within the other 8 approaches, estimate the percentage of these cDEGs identified by one approach - all seeds#####
library(tidyverse)
var.name <- c("BMI")
p.fun <- function(SMA.method.1 = NULL, SMA.method.2 = NULL, data.set = NULL) {
  d <- data.set |> 
    filter(.data[[SMA.method.1]] == 1) |> 
    summarise(n = n()) |> pull()
  n <- data.set |> 
    filter(.data[[SMA.method.1]] == 1, .data[[SMA.method.2]] == 1) |> 
    summarise(n = n()) |> pull()
  p <- n/d
  return(p)
}

SMA.method.all <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
p.final <- matrix(NA,
                  nrow = 8*9,
                  ncol = 3) |> data.frame()
colnames(p.final) <- c("approach_1", "approach_2", "percentage_cDEGs")

for (i in seq_along(SMA.method.all)) {
  SMA.method.vec <- SMA.method.all[-i]
  p.list <- vector(mode = "list", length = 10)
  names(p.list) <- seed.vec 
  for (z in seq_along(seed.vec)) {
    class.freq <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_%s.RDS", var.name, seed.vec[z]))
    class.freq <- class.freq[, -ncol(class.freq)]
    class.freq <- class.freq |> mutate(across(.cols = everything(), as.numeric))
    class.freq[, ncol(class.freq) + 1] <- rowSums(class.freq)
    colnames(class.freq)[ncol(class.freq)] <- "row.sum"
    p.vec <- vector(mode = "numeric", length = 8)
    for (k in seq_along(SMA.method.vec)) {
      p <- p.fun(SMA.method.1 = SMA.method.vec[k], SMA.method.2 = SMA.method.all[i], data.set = class.freq)
      p.vec[k] <- p
    }
    names(p.vec) <- SMA.method.vec
    p.list[[z]] <- p.vec
  }
  p.list.to.df <- do.call(rbind, p.list)
  p.list.to.df.avg <- p.list.to.df |> 
    as.data.frame() |> 
    summarise(across(.cols = everything(), mean, .names = "{.col}")) |>
    pivot_longer(cols = everything(), names_to = "approach_2", values_to = "percentage_cDEGs") |>
    mutate(approach_1 = SMA.method.all[i]) |> dplyr::select(approach_1, approach_2, percentage_cDEGs)
  p.final[(8*i-7):(8*i), ] <- p.list.to.df.avg
}


date.analysis <- format(Sys.Date(), "%Y%b%d")

theme_BMA <- function(
    base_size = 12,
    base_family = "Arial"
) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6), # Customize the box border in strip
      strip.placement = "inside"
    )
}

p.final <- mutate(p.final, percentage_cDEGs_2 = round(100*percentage_cDEGs))
p.final <- mutate(.data = p.final, across(.cols = starts_with("Approach"), .fns = ~ factor(.x, levels = SMA.method.all)))
g <- ggplot(data = p.final, 
       mapping = aes(x = approach_1, y = approach_2, color = percentage_cDEGs_2, size = percentage_cDEGs_2)) + 
  geom_point() + 
  theme_BMA() + 
  scale_color_viridis(direction = -1, option = "F") + 
  geom_text(mapping = aes(label = percentage_cDEGs_2), size = 4, nudge_x = 0.45, nudge_y = 0.25, color = "black") + 
  scale_size(range = c(1, 15)) + 
  labs(x = "Approach", y = "Approach", color = "Percentage \nof cDEGs") + 
  guides(size = "none") + 
  theme(legend.position = "right", 
        legend.key.height = unit(5, "line"),
        axis.text.y = element_text(color = "black"))

ggsave(filename = sprintf("../ApplicationResult/AddViz/shared_discovery/%s_%s_%s_dotplot.eps", date.analysis, var.name, "all_seed"),
       plot = g, device = cairo_ps, dpi = 600, width = 8, height = 8, units = "in")

p.final2 <- p.final |> filter(!(approach_1 %in% c("eBayes_MVM", "eBayes_UVM")), !(approach_2 %in% c("eBayes_MVM", "eBayes_UVM")))
g2 <- ggplot(data = p.final2, 
            mapping = aes(x = approach_1, y = approach_2, color = percentage_cDEGs_2, size = percentage_cDEGs_2)) + 
  geom_point() + 
  theme_BMA() + 
  scale_color_viridis(direction = -1, option = "F") + 
  geom_text(mapping = aes(label = percentage_cDEGs_2), size = 4, nudge_x = 0.45, nudge_y = 0.25, color = "black") + 
  scale_size(range = c(1, 15)) + # Also modify the color mapping to the mixmimum size
  labs(x = "Approach", y = "Approach", color = "Percentage \nof cDEGs") + 
  guides(size = "none") + 
  theme(legend.position = "right", 
        legend.key.height = unit(5, "line"),
        axis.text.y = element_text(color = "black"))
  
ggsave(filename = sprintf("../ApplicationResult/AddViz/shared_discovery/%s_%s_%s_dotplot2.eps", date.analysis, var.name, "all_seed"),
       plot = g2, device = cairo_ps, dpi = 600, width = 7, height = 7, units = "in")
