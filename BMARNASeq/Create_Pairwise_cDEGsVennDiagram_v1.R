easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "tidyverse", "edgeR", "DESeq2", "data.table", "ggVennDiagram", "gridExtra", "openxlsx", "cowplot") |> 
  suppressPackageStartupMessages()
Bayesfactor <- BMAseq:::Bayesfactor

threshold <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")

get_venndata <- function(seed.i = 8809678, 
                         var.name = "BMI", 
                         model.type.1 = "Multi", 
                         model.type.2 = "Uni",
                         model.method.1 = "BMAseq",
                         model.method.2 = "DESeq2",
                         threshold = 5000,
                         ...){
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
  return(list(cDEGs.1, cDEGs.2))
}


plot_venn <- function(seed.i = 8809678, 
                      var.name = "BMI", 
                      model.method.1 = "BMAseq",
                      model.method.2 = "DESeq2",
                      model.type.1 = "Multi",
                      model.type.2 = "Uni",
                      threshold = 5000,
                      ...) {

  pair_venndata <- get_venndata(seed.i = seed.i, 
                                var.name = var.name, 
                                model.method.1 = model.method.1,
                                model.method.2 = model.method.2,
                                model.type.1 = model.type.1,
                                model.type.2 = model.type.2,
                                threshold = threshold)
  
  # print(pair_venndata)
  
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
    model.type.2 == "Multi",
    paste0(model.method.2, "_MVM"),
    paste0(model.method.2, "_UVM")
  )))
  
  venn_data <- list(pair_venndata[[1]], pair_venndata[[2]])
  names(venn_data) <- label_names
  
  make_venn <- ggVennDiagram(
    x = venn_data,
    label = "count",
    label_alpha = 0,
    label_color = c("white", "white", "black")
  )
  make_venn$layers[[1]]$mapping <- aes(fill = name) # Set the discrete scale
  make_venn <- make_venn + 
    scale_color_manual(values = c("black", "black")) + # Black circle border
    scale_fill_manual(values = c("#3B499299", "#BB002199", "#00828099")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, colour = "red")) + 
  ggtitle(var.name)
  return(make_venn)
}

### One case
p1 <- plot_venn(seed.i = 8809678, 
                var.name = "BMI", 
                model.method.1 = "BMAseq",
                model.method.2 = "DESeq2",
                model.type.1 = "Multi",
                model.type.2 = "Uni",
                threshold = 2000)

p2 <- plot_venn(seed.i = 8809678, 
                var.name = "BMI", 
                model.method.1 = "BMAseq",
                model.method.2 = "DESeq2",
                model.type.1 = "Multi",
                model.type.2 = "Multi",
                threshold = 2000)

p3 <- plot_venn(seed.i = 8809678, 
                var.name = "BMI", 
                model.method.1 = "DESeq2",
                model.method.2 = "DESeq2",
                model.type.1 = "Multi",
                model.type.2 = "Uni",
                threshold = 2000)

plot_panel <- grid.arrange(p1, p2, p3, ncol = 3)
date.analysis <- format(Sys.Date(), "%Y%b%d")
file.name <- sprintf("../ApplicationResult/AddViz/VennDigram/Top5000/%s_%s.png", date.analysis, "BMI") # EPS does not support semi-transparency; tiff produces MB-size file
ggsave(filename = file.name, plot = plot_panel, device = "png", dpi = 600, width = 9, height = 6, units = "in")
