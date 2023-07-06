library(data.table)
library(tidyverse)
library(openxlsx)
var.vec <- c("BMI", "AGE", "SEX", "MHABNWBC")
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
threshold.vec <- 5000

# Approach 1
t0 <- microbenchmark::get_nanotime()
for (threshold.i in threshold.vec) {
  for (var.name in var.vec) {
    class.freq.list <- vector(mode = "list", length = length(seed.vec))
    names(class.freq.list) <- seed.vec
    for (seed.i in seed.vec) {
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/BMAseqMulti%s.RData", seed.i))
      names(BMAseq.eFDR.Main.train) = names(BMAseq.eFDR.Main.test) = var.vec
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/DESeq2Multi%s.RData", seed.i))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/edgeRMulti%s.RData", seed.i))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/eBayesMulti%s.RData", seed.i))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/MultiModel/VoomLimmaMulti%s.RData", seed.i))
      BMAseq.cDEGs <- intersect(names(BMAseq.eFDR.Main.train[[var.name]][1:threshold.i]), names(BMAseq.eFDR.Main.test[[var.name]][1:threshold.i]))
      DESeq2.MVM.cDEGs <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold.i], DESeq2.eFDR.GeneName.test[[var.name]][1:threshold.i])
      edgeR.MVM.cDEGs <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold.i], edgeR.eFDR.GeneName.test[[var.name]][1:threshold.i])
      eBayes.MVM.cDEGs <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold.i]), names(eBayes.eFDR.test2[[var.name]][1:threshold.i]))
      voom.limma.MVM.cDEGs <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold.i]), names(voom.eFDR.test2[[var.name]][1:threshold.i]))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/DESeq2Uni%s.RData", seed.i))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/edgeRUni%s.RData", seed.i))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/eBayesUni%s.RData", seed.i))
      load(sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/VoomLimmaUni%s.RData", seed.i))
      DESeq2.UVM.cDEGs <- intersect(DESeq2.eFDR.GeneName.train[[var.name]][1:threshold.i], DESeq2.eFDR.GeneName.test[[var.name]][1:threshold.i])
      edgeR.UVM.cDEGs <- intersect(edgeR.eFDR.GeneName.train[[var.name]][1:threshold.i], edgeR.eFDR.GeneName.test[[var.name]][1:threshold.i])
      eBayes.UVM.cDEGs <- intersect(names(eBayes.eFDR.train2[[var.name]][1:threshold.i]), names(eBayes.eFDR.test2[[var.name]][1:threshold.i]))
      voom.limma.UVM.cDEGs <- intersect(names(voom.eFDR.train2[[var.name]][1:threshold.i]), names(voom.eFDR.test2[[var.name]][1:threshold.i]))
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
t1 <- microbenchmark::get_nanotime()
t1 - t0
# [1] 2.775805e+12

                               
# Approach 2
t0 <- microbenchmark::get_nanotime()
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
t1 <- microbenchmark::get_nanotime()
t1 - t0
# [1] 2.974833e+12
