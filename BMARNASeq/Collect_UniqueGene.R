threshold <- 5000
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
date.analysis <- format(Sys.Date(), "%Y%b%d")
var.name <- "BMI"


cDEGs.type.vec <- c("BMAseq.cDEGs", 
                    "DESeq2.UVM.cDEGs",
                    "edgeR.UVM.cDEGs",
                    "eBayes.UVM.cDEGs",
                    "voom.limma.UVM.cDEGs",
                    "DESeq2.MVM.cDEGs",
                    "edgeR.MVM.cDEGs",
                    "eBayes.MVM.cDEGs",
                    "voom.limma.MVM.cDEGs")

for (var.name in var.vec) {
  for (target.type in cDEGs.type.vec) {
    unique.cDEGs.all.seed <- NULL
    target.type.idx <- which(cDEGs.type.vec == target.type)
    nontarget.type.idx <- which(cDEGs.type.vec != target.type)
    
    for (seed.i in seed.vec) {
      ##------BMAseq------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/Trim/BMAseqMulti%s.RData", seed.i) # Multivariable BMAseq model without interaction terms
      data.env <- local({load(file.name); environment()})
      names(data.env$BMAseq.eFDR.Main.train) = names(data.env$BMAseq.eFDR.Main.test) = var.vec
      BMAseq.cDEGs <- intersect(names(data.env$BMAseq.eFDR.Main.train[[var.name]][1:threshold]), 
                                names(data.env$BMAseq.eFDR.Main.test[[var.name]][1:threshold]))
      rm(list = "data.env")
      
      ##------DESeq2_UVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/Trim/DESeq2Uni%s.RData", seed.i) # Univariable DESeq2 model
      data.env <- local({load(file.name); environment()})
      DESeq2.UVM.cDEGs <- intersect(data.env$DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                                    data.env$DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
      rm(list = "data.env")
      
      ##------edgeR_UVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/Trim/edgeRUni%s.RData", seed.i) # Univariable edgeR model
      data.env <- local({load(file.name); environment()})
      edgeR.UVM.cDEGs <- intersect(data.env$edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                                   data.env$edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
      rm(list = "data.env")
      
      ##------eBayes_UVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/Trim/eBayesUni%s.RData", seed.i) # Univariable eBayes model
      data.env <- local({load(file.name); environment()})
      eBayes.UVM.cDEGs <- intersect(names(data.env$eBayes.eFDR.train2[[var.name]][1:threshold]), 
                                    names(data.env$eBayes.eFDR.test2[[var.name]][1:threshold]))
      rm(list = "data.env")
      
      ##------voom.limma_UVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/Trim/VoomLimmaUni%s.RData", seed.i) # Univariable VoomLimma model
      data.env <- local({load(file.name); environment()})
      voom.limma.UVM.cDEGs <- intersect(names(data.env$voom.eFDR.train2[[var.name]][1:threshold]), 
                                        names(data.env$voom.eFDR.test2[[var.name]][1:threshold]))
      rm(list = "data.env")
      
      ##------DESeq2_MVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/Trim/DESeq2Multi%s.RData", seed.i) # Multivariable DESeq2 model
      data.env <- local({load(file.name); environment()})
      DESeq2.MVM.cDEGs <- intersect(data.env$DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                                    data.env$DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
      rm(list = "data.env")
      
      ##------edgeR_MVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/Trim/edgeRMulti%s.RData", seed.i) # Multivariable edgeR model
      data.env <- local({load(file.name); environment()})
      edgeR.MVM.cDEGs <- intersect(data.env$edgeR.eFDR.GeneName.train[[var.name]][1:threshold],
                                   data.env$edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
      rm(list = "data.env")
      
      ##------eBayes_MVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/Trim/eBayesMulti%s.RData", seed.i) # Multivariable eBayes model
      data.env <- local({load(file.name); environment()})
      eBayes.MVM.cDEGs <- intersect(names(data.env$eBayes.eFDR.train2[[var.name]][1:threshold]), 
                                    names(data.env$eBayes.eFDR.test2[[var.name]][1:threshold]))
      rm(list = "data.env")
      
      ##------voom.limma_MVM------
      file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/Trim/VoomLimmaMulti%s.RData", seed.i) # Multivariable VoomLimma model
      data.env <- local({load(file.name); environment()})
      voom.limma.MVM.cDEGs <- intersect(names(data.env$voom.eFDR.train2[[var.name]][1:threshold]), 
                                        names(data.env$voom.eFDR.test2[[var.name]][1:threshold]))
      rm(list = "data.env")
      
      unique.cDEGs.per.seed <-
        setdiff(
          get(target.type),
          union(get(cDEGs.type.vec[nontarget.type.idx[1]]), get(cDEGs.type.vec[nontarget.type.idx[2]])) |>
            union(get(cDEGs.type.vec[nontarget.type.idx[3]])) |>
            union(get(cDEGs.type.vec[nontarget.type.idx[4]])) |>
            union(get(cDEGs.type.vec[nontarget.type.idx[5]])) |>
            union(get(cDEGs.type.vec[nontarget.type.idx[6]])) |>
            union(get(cDEGs.type.vec[nontarget.type.idx[7]])) |>
            union(get(cDEGs.type.vec[nontarget.type.idx[8]]))
        )
      
      method.name <- gsub("\\.", "_", 
                          gsub(".cDEGs", "", target.type))
      txt.file.name <- sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_uniq_cDEGs_%s_per_seed.txt", method.name, var.name)
      write(sprintf("\nThe following gene IDs are unique gene IDs associated with the main effect of %s identified by %s under the seed %s and threshold %s", 
                    var.name, method.name, seed.i, threshold), file = txt.file.name, append = T)
      write.table(unique.cDEGs.per.seed, 
                  file = txt.file.name, 
                  quote = F, col.names = F, append = T)
      
      unique.cDEGs.all.seed <- append(unique.cDEGs.all.seed, list(unique.cDEGs.per.seed))
    }
    names(unique.cDEGs.all.seed) <- seed.vec
    saveRDS(unique.cDEGs.all.seed, file = sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s.uniq.cDEGs.%s.all.seed.RDS", method.name, var.name))
  }
}
