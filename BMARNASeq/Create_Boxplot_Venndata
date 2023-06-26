easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", "openxlsx", "cowplot") |> 
  suppressPackageStartupMessages()
Bayesfactor <- BMAseq:::Bayesfactor

threshold.vec <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
model.method.vec <- c("DESeq2", "edgeR", "eBayes", "VoomLimma")

count_list_3 <- NULL
for (threshold.i in threshold.vec) {
  count_list_2 <- NULL
  for (var.i in var.vec) {
    count_list_1 <- NULL
    for (model.method.i in model.method.vec) {
      V1 = V2 = V3 = V4 = V5 = V6 = V7 = V8 = V9 = V10 = V11 = rep(NA, 7)
      count_df <- data.frame(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11)
      colnames(count_df) <- c("set.name", seed.vec)
      for (seed.i in seed.vec) { 
        tri_venndata <- get_venndata_call(seed.i = seed.i, 
                                     var.name = var.i, 
                                     model.method.1 = "BMAseq",
                                     model.method.2 = model.method.i,
                                     model.method.3 = model.method.i,
                                     model.type.1 = "Multi",
                                     model.type.2 = "Uni",
                                     model.type.3 = "Multi",
                                     threshold = threshold.i)
        model.method.1 = tri_venndata[[4]][["model.method.1"]]
        model.method.2 = model.method.i
        model.method.3 = model.method.i
        model.type.1 = tri_venndata[[4]][["model.type.1"]]
        model.type.2 = tri_venndata[[4]][["model.type.2"]]
        model.type.3 = tri_venndata[[4]][["model.type.3"]]
        
        label_names <- c(paste0(ifelse(
          model.method.1 == "BMAseq", "BMAseq",
          ifelse(model.type.1 == "Multi", 
                 paste0(model.method.1, "_MVM"),
                 paste0(model.method.1, "_UVM")
          )
        )),
        paste0(ifelse(
          model.method.2 == "BMAseq", "BMAseq",
          ifelse(model.type.2 == "Multi", 
                 paste0(model.method.2, "_MVM"),
                 paste0(model.method.2, "_UVM")
          )
        )),
        paste0(ifelse(
          model.method.3 == "BMAseq", "BMAseq",
          ifelse(model.type.3 == "Multi",
                 paste0(model.method.3, "_MVM"),
                 paste0(model.method.3, "_UVM")
          )))
        )
        label_names <- ifelse(grepl("VoomLimma", label_names) == T, 
                              gsub("VoomLimma", "voom.limma", label_names), 
                              label_names)
        
        tri_venndata <- list(tri_venndata[[1]], tri_venndata[[2]], tri_venndata[[3]])
        names(tri_venndata) <- label_names
        processed_data <- Venn(tri_venndata) |> process_data()
        count_data <- processed_data@region |> as.data.frame() |> select(name, count)
        col_idx <- Position(function(x) x == seed.i, seed.vec) # Locate the position of a seed in a seed vector
        if (any(is.na(count_df[, 1]))) { # Assign the set name only once to save the computational steps
          count_df[, 1] <- count_data[["name"]]
        }
        count_df[, col_idx + 1] <- count_data[["count"]]
        # merge(count_df, count_data, by = "name")
      }
      count_list_per_method <- list(count_df)
      names(count_list_per_method) <- paste(threshold.i, var.i, model.method.i, sep = "+")
      count_list_1 <- c(count_list_1, count_list_per_method)
    }
    count_list_per_var <- list(count_list_1)
    names(count_list_per_var) <- var.i
    count_list_2 <- c(count_list_2, count_list_per_var)
  }
  count_list_per_threshold <- list(count_list_2)
  names(count_list_per_threshold) <- threshold.i
  count_list_3 <- c(count_list_3, count_list_per_threshold)
  return(count_list_3)
}

analysis.date <- format(Sys.Date(), "%Y%b%d")
saveRDS(object = count_list_3, file = paste(analysis.date, "count_list_3.RDS", sep = "_"))

##------Make a boxplot data organizing function------
get_boxplotdata <- function(var.name = var.i) {
  df <- do.call(rbind, count_list_3[["5000"]][[var.name]])
  df$model.method <- c(rep(c("DESeq2", "edgeR", "eBayes", "voom.limma"), each = 7))
  df.long <- pivot_longer(data = df,
                          cols = !c(set.name, model.method),
                          names_to = "seed",
                          values_to = "count")
  df.long$set.name <- ifelse(df.long$set.name == "BMAseq", 
                             paste0(df.long$set.name, "_vs_", df.long$model.method),
                             df.long$set.name)
  df.long$set.name <- gsub("^BMAseq_vs_DESeq2$", "BMAseq_vs_DESeq2_UVM_vs_DESeq2_MVM", df.long$set.name)
  df.long$set.name <- gsub("^BMAseq_vs_edgeR$", "BMAseq_vs_edgeR_UVM_vs_edgeR_MVM", df.long$set.name)
  df.long$set.name <- gsub("^BMAseq_vs_eBayes$", "BMAseq_vs_eBayes_UVM_vs_eBayes_MVM", df.long$set.name)
  df.long$set.name <- gsub("^BMAseq_vs_voom.limma$", "BMAseq_vs_voom.limma_UVM_vs_voom.limma_MVM", df.long$set.name)
  df.long$set.name <- gsub("^DESeq2_UVM$", "DESeq2_UVM_vs_BMAseq_vs_DESeq2_MVM", df.long$set.name)
  df.long$set.name <- gsub("^DESeq2_MVM$", "DESeq2_MVM_vs_BMAseq_vs_DESeq2_UVM", df.long$set.name)
  df.long$set.name <- gsub("^edgeR_UVM$", "edgeR_UVM_vs_BMAseq_vs_edgeR_MVM", df.long$set.name)
  df.long$set.name <- gsub("^edgeR_MVM$", "edgeR_MVM_vs_BMAseq_vs_edgeR_UVM", df.long$set.name)
  df.long$set.name <- gsub("^eBayes_UVM$", "eBayes_UVM_vs_BMAseq_vs_eBayes_MVM", df.long$set.name)
  df.long$set.name <- gsub("^eBayes_MVM$", "eBayes_MVM_vs_BMAseq_vs_eBayes_UVM", df.long$set.name)
  df.long$set.name <- gsub("^voom.limma_UVM$", "voom.limma_UVM_vs_BMAseq_vs_voom.limma_MVM", df.long$set.name)
  df.long$set.name <- gsub("^voom.limma_MVM$", "voom.limma_MVM_vs_BMAseq_vs_voom.limma_UVM", df.long$set.name)
  
  df.long$set.name <- factor(df.long$set.name, 
                             levels = c(
                               "BMAseq..DESeq2_UVM..DESeq2_MVM",
                               "BMAseq..edgeR_UVM..edgeR_MVM",
                               "BMAseq..eBayes_UVM..eBayes_MVM",
                               "BMAseq..voom.limma_UVM..voom.limma_MVM",
                               
                               "BMAseq..DESeq2_UVM",
                               "BMAseq..DESeq2_MVM",
                               "BMAseq..edgeR_UVM",
                               "BMAseq..edgeR_MVM",
                               "BMAseq..eBayes_UVM",
                               "BMAseq..eBayes_MVM",
                               "BMAseq..voom.limma_UVM",
                               "BMAseq..voom.limma_MVM",
                               "DESeq2_UVM..DESeq2_MVM",
                               "edgeR_UVM..edgeR_MVM",
                               "eBayes_UVM..eBayes_MVM",
                               "voom.limma_UVM..voom.limma_MVM",
                               
                               "BMAseq_vs_DESeq2_UVM_vs_DESeq2_MVM",
                               "BMAseq_vs_edgeR_UVM_vs_edgeR_MVM",
                               "BMAseq_vs_eBayes_UVM_vs_eBayes_MVM",
                               "BMAseq_vs_voom.limma_UVM_vs_voom.limma_MVM",
                               "DESeq2_UVM_vs_BMAseq_vs_DESeq2_MVM",
                               "DESeq2_MVM_vs_BMAseq_vs_DESeq2_UVM",
                               "edgeR_UVM_vs_BMAseq_vs_edgeR_MVM",
                               "edgeR_MVM_vs_BMAseq_vs_edgeR_UVM",
                               "eBayes_UVM_vs_BMAseq_vs_eBayes_MVM",
                               "eBayes_MVM_vs_BMAseq_vs_eBayes_UVM",
                               "voom.limma_UVM_vs_BMAseq_vs_voom.limma_MVM",
                               "voom.limma_MVM_vs_BMAseq_vs_voom.limma_UVM"
                             ))
  return(df.long)
}

##------Set the plot theme------
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

##------Make a boxplot plotting function------
make_boxplot <- function(var.name = var.name, data.set = data.set) {
  p <- ggplot(data = data.set,
              aes(y = count, x = set.name)) + 
    geom_boxplot(aes(group = set.name, color = set.name),
                 outlier.shape = NA) + 
    geom_jitter(aes(color = set.name), size = 0.8) + 
    geom_vline(xintercept = 16.5, linetype = "dashed") + 
    annotate(geom = "text", x = 8.5, y = 1200, label = "Shared discovery", size = 5) + 
    annotate(geom = "text", x = 22.5, y = 1200, label = "Unique discovery", size = 5) + 
    xlab("Set name") + 
    ylab(paste0("Number of shared ", var.name, "-related cDEGs")) + 
    theme_BMA() + 
    scale_color_manual(values = (c(rep("#DF8F44FF", 4), rep("#B24745FF", 12), rep("#6A6599FF", 12)))) + 
    coord_flip() 
  date.analysis <- format(Sys.Date(), "%Y%b%d")
  ggsave(filename = sprintf("../ApplicationResult/AddViz/VennDigram/Top5000/All/%s_%s.png", date.analysis, var.name),
         plot = p, device = "png", dpi = 600, width = 11, height = 7, units = "in")
}

var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
for (var.i in var.vec) {
  boxplot.data <- get_boxplotdata(var.name = var.i)
  make_boxplot(var.name = var.i, data.set = boxplot.data)
}
