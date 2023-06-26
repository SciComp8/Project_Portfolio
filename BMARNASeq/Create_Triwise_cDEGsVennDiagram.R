easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", "openxlsx", "cowplot") |> 
  suppressPackageStartupMessages()
Bayesfactor <- BMAseq:::Bayesfactor

threshold.vec <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
model.method.vec <- c("DESeq2", "edgeR", "eBayes", "VoomLimma")

get_venndata <- function(seed.i = 8809678, 
                         var.name = "BMI", 
                         model.type.1 = "Multi", 
                         model.type.2 = "Uni",
                         model.type.3 = "Multi",
                         model.method.1 = "BMAseq",
                         model.method.2 = "DESeq2",
                         model.method.3 = "DESeq2",
                         threshold = 5000,
                         ...) {
  file.name.1 <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                         model.type.1, model.method.1, model.type.1, seed.i)
  load(file.name.1)
  if (model.method.1 == "BMAseq") {
    names(BMAseq.eFDR.Main.train) = names(BMAseq.eFDR.Main.test) = var.vec
    cDEGs.1 <- intersect(names(BMAseq.eFDR.Main.train[[var.name]][1:threshold]), 
                         names(BMAseq.eFDR.Main.test[[var.name]][1:threshold]))
  } else if (model.method.1 == "DESeq2") {
    cDEGs.1 <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                         DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.1 == "edgeR") {
    cDEGs.1 <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                         edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.1 == "eBayes") {
    cDEGs.1 <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold]),
                         names(eBayes.eFDR.test2[[var.name]][1:threshold]))
  } else if (model.method.1 == "VoomLimma") {
    cDEGs.1 <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold]),
                         names(voom.eFDR.test2[[var.name]][1:threshold]))
  }
  
  file.name.2 <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                         model.type.2, model.method.2, model.type.2, seed.i)
  load(file.name.2)
  if (model.method.2 == "DESeq2") {
    cDEGs.2 <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                         DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.2 == "edgeR") {
    cDEGs.2 <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                         edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.2 == "eBayes") {
    cDEGs.2 <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold]),
                         names(eBayes.eFDR.test2[[var.name]][1:threshold]))
  } else if (model.method.2 == "VoomLimma") {
    cDEGs.2 <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold]),
                         names(voom.eFDR.test2[[var.name]][1:threshold]))
  }
  
  file.name.3 <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                         model.type.3, model.method.3, model.type.3, seed.i)
  load(file.name.3)
  if (model.method.3 == "DESeq2") {
    cDEGs.3 <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                         DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.3 == "edgeR") {
    cDEGs.3 <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                         edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
  } else if (model.method.3 == "eBayes") {
    cDEGs.3 <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold]),
                         names(eBayes.eFDR.test2[[var.name]][1:threshold]))
  } else if (model.method.3 == "VoomLimma") {
    cDEGs.3 <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold]),
                         names(voom.eFDR.test2[[var.name]][1:threshold]))
  }
  
  return(list(cDEGs.1, cDEGs.2, cDEGs.3))
}


plot_venn <- function(seed.i = 8809678, 
                      var.name = "BMI", 
                      model.method.1 = "BMAseq",
                      model.method.2 = "DESeq2",
                      model.method.3 = "DESeq2",
                      model.type.1 = "Multi",
                      model.type.2 = "Uni",
                      model.type.3 = "Multi",
                      threshold = 5000,
                      ...) {

  tri_venndata <- get_venndata(seed.i = seed.i, 
                               var.name = var.name, 
                               model.method.1 = model.method.1,
                               model.method.2 = model.method.2,
                               model.method.3 = model.method.3,
                               model.type.1 = model.type.1,
                               model.type.2 = model.type.2,
                               model.type.3 = model.type.3,
                               threshold = threshold)
  
  # print(tri_venndata)
  
  label_names <- c(paste0(ifelse(
    model.method.1 == "BMAseq",
    "BMAseq",
    ifelse(
      model.type.1 == "Multi",
      paste0(model.method.1, "_MVM"),
      paste0(model.method.1, "_UVM")
    )
  )),
  paste0(ifelse(
    model.method.2 == "BMAseq",
    "BMAseq",
    ifelse(
      model.type.2 == "Multi",
      paste0(model.method.2, "_MVM"),
      paste0(model.method.2, "_UVM")
    )
  )),
  paste0(ifelse(
    model.method.3 == "BMAseq",
    "BMAseq",
    ifelse(
      model.type.3 == "Multi",
      paste0(model.method.3, "_MVM"),
      paste0(model.method.3, "_UVM")
    )
  )))
  label_names <- ifelse(grepl("VoomLimma", label_names) == T, gsub("VoomLimma", "voom.limma", label_names), label_names)
  
  venn_data <- list(tri_venndata[[1]], tri_venndata[[2]], tri_venndata[[3]])
  names(venn_data) <- label_names
  
  make_venn <- ggVennDiagram(
    x = venn_data,
    label = "count",
    label_alpha = 0,
    label_color = "black",
    label_size = 10,
    set_size = 6
  )
  make_venn$layers[[1]]$mapping <- aes(fill = name) # Set the discrete scale
  make_venn <- make_venn + 
    scale_x_continuous(expand = expansion(mult = 0.6)) + # Expand axis to show long set labels
    scale_color_manual(values = c("black", "black", "black")) + # Black circle border
    scale_fill_manual(values = c("#3B499299", "#BB002199", "#00828099", "#A2005699", "#80818099", "#63197999", "#008B4599")) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5, colour = "red")) + 
  ggtitle(paste0(var.name, " (seed ", seed.i, ", ranking threshold ", threshold, ")"))
  return(make_venn)
}

# Loop every seed, every var.name, every ranking threshold
for (threshold.i in threshold.vec) {
  for (var.i in var.vec) {
    for (seed.i in seed.vec) {
      for (model.method.i in model.method.vec) {
        p1 <- plot_venn(seed.i = seed.i, 
                        var.name = var.i, 
                        model.method.1 = "BMAseq",
                        model.method.2 = model.method.i,
                        model.method.3 = model.method.i,
                        model.type.1 = "Multi",
                        model.type.2 = "Uni",
                        model.type.3 = "Multi",
                        threshold = threshold.i)
        date.analysis <- format(Sys.Date(), "%Y%b%d")
        file.name <- sprintf("../ApplicationResult/AddViz/VennDigram/Top5000/%s/%s_%s_%s.png", model.method.i, date.analysis, var.i, seed.i) 
        ggsave(filename = file.name, plot = p1, device = "png", dpi = 600, width = 12, height = 9, units = "in")        
      }
    }
  }
}
