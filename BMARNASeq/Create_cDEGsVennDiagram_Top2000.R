library(ggVennDiagram)
threshold <- 2000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE","SEX", "MHABNWBC", "BMIxSEX")

date.analysis <- format(Sys.Date(), "%Y%b%d")
make_vennplot <-function(var.name = NULL, threshold = NULL) {
  if (var.name == "BMIxSEX") {
    var.name <- "BMI"
    BMAseq.cDEGs <- intersect(names(BMAseq.eFDR.Interaction.train[[var.name]][1:threshold]), names(BMAseq.eFDR.Interaction.test[[var.name]][1:threshold]))
  } else {
    BMAseq.cDEGs <- intersect(names(BMAseq.eFDR.Main.train[[var.name]][1:threshold]), names(BMAseq.eFDR.Main.test[[var.name]][1:threshold]))
  }
  DESeq2.cDEGs <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
  edgeR.cDEGs <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold], edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
  eBayes.cDEGs <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold]), names(eBayes.eFDR.test2[[var.name]][1:threshold]))
  voom.limma.cDEGs <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold]), names(voom.eFDR.test2[[var.name]][1:threshold]))
  g <- ggVennDiagram(list(BMAseq = BMAseq.cDEGs, 
                          DESeq2 = DESeq2.cDEGs,
                          edgeR = edgeR.cDEGs,
                          eBayes = eBayes.cDEGs,
                          voom.limma = voom.limma.cDEGs),
                     label_alpha = 0, label_color = "white") +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, colour = "red", size = 18)) + 
    ggtitle(paste("The intersection of common differentially expressed genes \nassociated with the", var.name)) + 
    scale_fill_distiller(palette = "Dark2", direction = -1)
  
  ggsave(filename = paste0("../ApplicationResult/Multi_Interaction/RandomSeed/TMM_Top2000/Vennplot/", date.analysis, "_", var.name, "_", seed.i, ".png"),
         plot = g,
         device = "png",
         width = 10,
         height = 7,
         units = "in")
}

for (seed.i in seed.vec) {
  load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/BMAseqMultiInt%s.RData", seed.i))
  load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/DESeq2MultiInt%s.RData", seed.i))
  load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/edgeRMultiInt%s.RData", seed.i))
  load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/eBayesMultiInt%s.RData", seed.i))
  load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/VoomLimmaMultiInt%s.RData", seed.i))
  for (var.i in var.vec) {
    make_vennplot(var.name = var.i, threshold = 2000)
  }
}
