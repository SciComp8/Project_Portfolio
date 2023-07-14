library(clusterProfiler)
library(org.Hs.eg.db)

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

ora.obj.list <- vector(mode = "list", length = 9)
names(ora.obj.list) <- method.type.vec
ora.obj.list <- list(ora.obj.list, ora.obj.list, ora.obj.list, ora.obj.list)
names(ora.obj.list) <- var.vec

for (var.name in var.vec) {
  for (method.name in method.type.vec) {
    if (grepl("voom.limma", method.name)) {
      method.name <- paste0("voom_limma", sub(".*_", "_", method.name))
    }
    unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s.uniq.cDEGs.%s.all.seed.RDS",
                                             method.name, var.name))
    unique.cDEGs <- unique.cDEGs.all.seed[[as.character(seed.i)]]
    unique.cDEGs.8809678 <- sub("\\..*", "", unique.cDEGs) # Delete any character after the comma
    
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
    ora.obj.list[[var.name]][[method.name]] <- ora.obj.ordered
  }
}
