t0 <- microbenchmark::get_nanotime()
cDEGs.all <- Reduce(union, list(BMAseq.cDEGs, DESeq2.MVM.cDEGs, edgeR.MVM.cDEGs, eBayes.MVM.cDEGs, voom.limma.MVM.cDEGs,
                                DESeq2.UVM.cDEGs, edgeR.UVM.cDEGs, eBayes.UVM.cDEGs, voom.limma.UVM.cDEGs)) 
t1 <- microbenchmark::get_nanotime()
t1-t0

# 2x faster
t0 <- microbenchmark::get_nanotime()
cDEGs.all.2 <- unique(c(BMAseq.cDEGs, DESeq2.MVM.cDEGs, edgeR.MVM.cDEGs, eBayes.MVM.cDEGs, voom.limma.MVM.cDEGs,
                        DESeq2.UVM.cDEGs, edgeR.UVM.cDEGs, eBayes.UVM.cDEGs, voom.limma.UVM.cDEGs))
t1 <- microbenchmark::get_nanotime()
t1-t0
