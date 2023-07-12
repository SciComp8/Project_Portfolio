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
      data.set <- mutate(data.set, qscore = -log(p.adjust, base = 10))
      color.upper <- range(data.set$qscore)[2] |> round(2)
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
         color = ifelse(color.type == "pvalue", "p-value", "q score")) +
    ggtitle(paste0("Unique cDEGs of ", var.name)) + 
    theme_bw(base_size = 14, base_family = "Arial") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 13, color = "black")) + 
    guides(size = guide_legend(order = 1)) # Order the discrete legend
  
  return(p)
}


date.analysis <- format(Sys.Date(), "%Y%b%d")
seed.i <- 8809678
threshold <- 5000
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC") 
method.type.vec <- c("BMAseq", 
                     "DESeq2_UVM",
                     "edgeR_UVM",
                     "eBayes_UVM",
                     "voom.limma_UVM",
                     "DESeq2_MVM",
                     "edgeR_MVM",
                     "eBayes_MVM",
                     "voom.limma_MVM")


for (var.name in var.vec) {
  for (method.name in method.type.vec) {
    if (grepl("voom.limma", method.name)) {
      method.name <- paste0("voom_limma", sub(".*_", "_", method.name))
    }
    unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s.uniq.cDEGs.%s.all.seed.RDS",
                                             method.name, var.name))
    unique.cDEGs <- unique.cDEGs.all.seed[[as.character(seed.i)]]
    unique.cDEGs.8809678 <- sub("\\..*", "", unique.cDEGs)      # Delete any character after the comma

    ora.obj <-
      enrichGO(
        gene          = unique.cDEGs.8809678,
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
    
    p <- make_dotplot(data.set = ora.obj.ordered, 
                      top.n = 20, 
                      method.name = method.i, 
                      color.lower = NULL, 
                      color.upper = NULL, 
                      color.type = "qscore")
    
    ggsave(filename = sprintf("../ApplicationResult/AddViz/DotPlot/%s_%s_%s_%s.eps", date.analysis, method.name, var.name, seed.i),
           plot = p, device = cairo_ps, dpi = 600, width = 10, height = 6, units = "in")
    
  }
}
