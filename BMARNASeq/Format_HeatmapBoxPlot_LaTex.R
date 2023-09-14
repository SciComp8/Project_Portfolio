##------Set up the environment------
library(data.table)
library(tidyverse)
library(openxlsx)
library(viridis)
library(grDevices)
source("theme_BMA.R")
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.vec <- "BMI"
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- 5000
threshold.i <- 5000

# Make a function to count the number of 1's in a string
count_ones <- function(string) {
  sum(as.integer(strsplit(string, "")[[1]]))
}

# Make a function to calculate the number of ones and their positions in a string
get_ones_info <- function(str) {
  ones_count <- sum(strsplit(str, "")[[1]] == "1")
  ones_positions <- which(strsplit(str, "")[[1]] == "1")
  return(list(string = str, count = ones_count, positions = ones_positions))
}

##------Estimate the number of cDEGs of each class per seed------
for (threshold.i in threshold.vec) {
  for (var.name in var.vec) {
    class.freq.list <- vector(mode = "list", length = length(seed.vec))
    names(class.freq.list) <- seed.vec
    for (seed.i in seed.vec) {
      target.files <- list.files(path = "../ApplicationData/derived/RandomSeed/Top5000/MultiModel/", pattern = sprintf("*Multi%s.RData", seed.i), full.names = T)
      temp <- new.env()
      lapply(target.files, load, temp)
      names(temp$BMAseq.eFDR.Main.train) = names(temp$BMAseq.eFDR.Main.test) = var.vec
      BMAseq.cDEGs <- with(temp, intersect(names(BMAseq.eFDR.Main.train[[var.name]][1:threshold.i]), names(BMAseq.eFDR.Main.test[[var.name]][1:threshold.i])))
      DESeq2.MVM.cDEGs <- with(temp, intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold.i], DESeq2.eFDR.GeneName.test[[var.name]][1:threshold.i]))
      edgeR.MVM.cDEGs <- with(temp, intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold.i], edgeR.eFDR.GeneName.test[[var.name]][1:threshold.i]))
      eBayes.MVM.cDEGs <- with(temp, intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold.i]), names(eBayes.eFDR.test2[[var.name]][1:threshold.i])))
      voom.limma.MVM.cDEGs <- with(temp, intersect(names(voom.eFDR.train2[[var.name]][1:threshold.i]), names(voom.eFDR.test2[[var.name]][1:threshold.i])))
      rm(temp)
      
      target.files <- list.files(path = "../ApplicationData/derived/RandomSeed/Top5000/UniModel/", pattern = sprintf("*Uni%s.RData", seed.i), full.names = T)
      temp <- new.env()
      lapply(target.files, load, temp)
      DESeq2.UVM.cDEGs <- with(temp, intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold.i], DESeq2.eFDR.GeneName.test[[var.name]][1:threshold.i]))
      edgeR.UVM.cDEGs <- with(temp, intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold.i], edgeR.eFDR.GeneName.test[[var.name]][1:threshold.i]))
      eBayes.UVM.cDEGs <- with(temp, intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold.i]), names(eBayes.eFDR.test2[[var.name]][1:threshold.i])))
      voom.limma.UVM.cDEGs <- with(temp, intersect(names(voom.eFDR.train2[[var.name]][1:threshold.i]), names(voom.eFDR.test2[[var.name]][1:threshold.i])))
      rm(temp)
      cDEGs.all <- unique(c(BMAseq.cDEGs, DESeq2.MVM.cDEGs, edgeR.MVM.cDEGs, eBayes.MVM.cDEGs, voom.limma.MVM.cDEGs,
                            DESeq2.UVM.cDEGs, edgeR.UVM.cDEGs, eBayes.UVM.cDEGs, voom.limma.UVM.cDEGs))
      no.cDEGs <- length(cDEGs.all)
      cDEGs.mat <- matrix(0, nrow = no.cDEGs, ncol = 10) 
      rownames(cDEGs.mat) <- cDEGs.all
      colnames(cDEGs.mat) <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM",
                               "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM", "Class")
      cDEGs.mat[, 1] <- as.integer(cDEGs.all %in% BMAseq.cDEGs)
      cDEGs.mat[, 2] <- as.integer(cDEGs.all %in% DESeq2.UVM.cDEGs)
      cDEGs.mat[, 3] <- as.integer(cDEGs.all %in% DESeq2.MVM.cDEGs)
      cDEGs.mat[, 4] <- as.integer(cDEGs.all %in% edgeR.UVM.cDEGs)
      cDEGs.mat[, 5] <- as.integer(cDEGs.all %in% edgeR.MVM.cDEGs)
      cDEGs.mat[, 6] <- as.integer(cDEGs.all %in% eBayes.UVM.cDEGs)
      cDEGs.mat[, 7] <- as.integer(cDEGs.all %in% eBayes.MVM.cDEGs)
      cDEGs.mat[, 8] <- as.integer(cDEGs.all %in% voom.limma.UVM.cDEGs)
      cDEGs.mat[, 9] <- as.integer(cDEGs.all %in% voom.limma.MVM.cDEGs)
      cDEGs.mat[, 10] <- apply(cDEGs.mat[, 1:9], 1, function(x) paste(x, collapse = ""))
      cDEGs.df <- cDEGs.mat |> as.data.frame()
      
      file.name <- paste0("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/", var.name, "_", threshold.i, "_", seed.i)
      saveRDS(cDEGs.df, paste0(file.name, ".RDS"))
      write.xlsx(x = cDEGs.df, file = paste0(file.name, ".xlsx"), rowNames = T, colWidths = 40, firstActiveRow = 2, firstActiveCol = 2)
      
      class.freq <- cDEGs.df |>
        dplyr::count(Class, sort = T, name = "Frequency") |>
        mutate(Seed = seed.i) |>
        select(Seed, Class, Frequency)
      
      class.freq.list[[seed.i]] <- class.freq
    }
    class.freq.df <- do.call(rbind, class.freq.list)
    saveRDS(class.freq.df, file = paste0("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/", var.name, "_", threshold.i, "_ClassFreq"))
  }
}

##------Make the integrative heatmap boxplot------
for (var.name in var.vec) {
  class.freq <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_ClassFreq", var.name))
  class.freq.top10.per.seed <- class.freq |> 
    group_by(Seed) |> 
    slice_max(order_by = Frequency, n = 10, with_ties = F) |>
    ungroup() 
  
  # Filter out the top 10 most frequent classes that appear less than 5 times across 10 seeds
  uniq_class <- unique(class.freq.top10.per.seed$Class)
  uniq_class_freq <- table(factor(class.freq.top10.per.seed$Class, levels = uniq_class))
  uniq_class_eligible <- uniq_class[uniq_class_freq >= 5]
  class.freq.top10.per.seed <- class.freq.top10.per.seed |> filter(Class %in% uniq_class_eligible)
  
  # Factorize and order the levels based on count of 1's
  uniq_class_2 <- unique(class.freq.top10.per.seed$Class)
  
  ones_info <- lapply(uniq_class_2, get_ones_info)
  
  df <- do.call(rbind, ones_info) |> data.frame()
  df <- df[order(df$count |> unlist(), sapply(df$positions, function(x) min(x))), ]
  new_class_level <- df$string |> unlist()
  
  # Reorder the class levels for plotting
  class.freq.top10.per.seed$Class <- factor(class.freq.top10.per.seed$Class, levels = rev(new_class_level))
  
  p1 <- ggplot(data = class.freq.top10.per.seed,
               aes(y = Frequency,
                   x = Class)) + 
    geom_boxplot(aes(group = Class, color = Class),
                 outlier.shape = NA) + 
    geom_jitter(aes(color = Class), size = 0.8) + 
    ylab(paste0("Number of cDEGs")) + 
    xlab("") + 
    theme_BMA() + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 0),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    scale_color_viridis(discrete = T,
                        direction = -1)
  
  aux <- ggplot_build(p1)
  y_ori <- aux[["layout"]][["panel_scales_y"]][[1]][["range"]][["range"]][2]
  if (var.name == "BMI") {
    y_pos <- y_ori - 400 # Set the position of label along the x axis
  } else if (var.name == "AGE") {
    y_pos <- y_ori - 400
  } else if (var.name == "SEX") {
    y_pos <- y_ori - 400
  } else if (var.name == "MHABNWBC") {
    y_pos <- y_ori - 900
  }
  
  p1 <- p1 + 
    annotate("text", x = 9, y = y_pos, label = "m1: BMAseq", hjust = 0) +
    annotate("text", x = 8.5, y = y_pos, label = "m2: DESeq2_UVM", hjust = 0) +
    annotate("text", x = 8, y = y_pos, label = "m3: DESeq2_MVM", hjust = 0) +
    annotate("text", x = 7.5, y = y_pos, label = "m4: edgeR_UVM", hjust = 0) +
    annotate("text", x = 7, y = y_pos, label = "m5: edgeR_MVM", hjust = 0) +
    annotate("text", x = 6.5, y = y_pos, label = "m6: eBayes_UVM", hjust = 0) +
    annotate("text", x = 6, y = y_pos, label = "m7: eBayes_MVM", hjust = 0) +
    annotate("text", x = 5.5, y = y_pos, label = "m8: voom.limma_UVM", hjust = 0) +
    annotate("text", x = 5, y = y_pos, label = "m9: voom.limma_MVM", hjust = 0)
  
  # Set the input vector
  input_vector <- uniq_class_2
  
  # Set the approach vector
  approach_positions <- c("BMAseq", "DESeq2_UVM", "DESeq2_MVM", "edgeR_UVM", "edgeR_MVM", "eBayes_UVM", "eBayes_MVM", "voom.limma_UVM", "voom.limma_MVM")
  
  # Set the empty data frame to store the results
  result_df <- data.frame(Y = character(), X = character(), Z = integer(), stringsAsFactors = F)
  
  # Process each string in the input vector
  for (i in seq_along(input_vector)) {
    string <- input_vector[i]
    
    # Find the positions of 1s in the string
    ones_positions <- which(strsplit(string, "")[[1]] == "1")
    
    # Create rows for the data frame based on the positions
    for (pos in ones_positions) {
      # Map the position of 1 to the approach
      row <- data.frame(Y = string, X = approach_positions[pos], Z = 1)
      result_df <- rbind(result_df, row)
    }
  }
  
  result_df$X <- factor(result_df$X, levels = approach_positions)
  result_df$Y <- factor(result_df$Y, levels = rev(new_class_level))
  p2 <- ggplot(result_df, aes(X, Y, fill = Z)) + 
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1,
              aes(fill = Y)) + 
    labs(x = "Method",
         y = "Detection way") + 
    theme_BMA() + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1, color = "black", vjust = 0.5),
          axis.text.y = element_text(color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)) + 
    scale_x_discrete(labels = paste0("m", 1:9)) + # Relabel the x ticks
    scale_fill_viridis(discrete = T,
                       direction = -1)
  
  g <- cowplot::plot_grid(p2, p1)
  ggsave(filename = sprintf("../ApplicationResult/AddViz/HeatmapBoxplot/%s_%s_%s_latex.eps", date.analysis, var.name, threshold.i),
         plot = g,
         device = cairo_ps,
         dpi = 600,
         width = 8,
         height = 4,
         units = "in")
}
