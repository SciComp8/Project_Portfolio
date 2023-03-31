 method <- "VoomLimma"
 if (!(method %in% c("BMAseq", "DESeq2", "edgeR"))){
    stop("The method is unavailable. It should be either one of BMAseq, DESeq2, and edgeR.")
  }
