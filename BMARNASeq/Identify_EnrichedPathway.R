library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
keytypes(org.Hs.eg.db)

##------Load image------
save.object <- ls(pattern = "*.go.ora") |> mget() # Return a list
seed.num = 8809678
date.analysis <- format(Sys.Date(), "%Y%b%d")
save(save.object, 
     file = sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_EnrichResult_%s.RData", date.analysis, seed.num))
file.name <- "../ApplicationResult/UniqueGene/TMM_Top5000/2023Jun30_EnrichResult_8809678.RData"
data.env <- local({load(file.name); environment()}) # Load image in a different environment

##------Load cDEGs data------
seed.i = 8809678
var.name = "BMI"
threshold = 5000

##------BMAseq-----
method.name = "BMAseq"
model.type = "Multi"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
BMAseq.cDEGs <- intersect(names(data.env$BMAseq.eFDR.Main.train[[var.name]][1:threshold]), 
                          names(data.env$BMAseq.eFDR.Main.test[[var.name]][1:threshold]))
BMAseq.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s.uniq.cDEGs.%s.all.seed.RDS",
                                                method.name, var.name))
BMAseq.unique.cDEGs <- BMAseq.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
BMAseq.cDEGs.8809678 <- sub("\\..*", "", BMAseq.cDEGs) 
BMAseq.unique.cDEGs.8809678 <- sub("\\..*", "", BMAseq.unique.cDEGs) 

# GO over-representation analysis
BMAseq.go.ora <-
  enrichGO(
    gene          = BMAseq.unique.cDEGs.8809678,
    # universe      = BMAseq.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
BMAseq.go.ora@result |> View()  
BMAseq.go.ora@result[order(BMAseq.go.ora@result$Count, decreasing = T), ] |> View()
BMAseq.go.ora@result[order(BMAseq.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "small molecule catabolic process"               
# [2] "amino acid metabolic process"                   
# [3] "organic acid catabolic process"                 
# [4] "carboxylic acid catabolic process"              
# [5] "axonogenesis"                                   
# [6] "purine nucleotide biosynthetic process"         
# [7] "purine-containing compound biosynthetic process"
# [8] "nucleotide biosynthetic process"                
# [9] "nucleoside phosphate biosynthetic process"      
# [10] "lipid catabolic process"          
BMAseq.go.ora@result[order(BMAseq.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "MTHFS/APOBEC3H/ECHS1/RIDA/PPARA/APOBEC3G/DBT/ALDH5A1/AUH/DECR1/GCSH"
# [2] "MTHFS/ECHS1/RIDA/NARS1/MARS2/DBT/ALDH5A1/APIP/AUH/GCSH"             
# [3] "MTHFS/ECHS1/RIDA/PPARA/DBT/ALDH5A1/AUH/DECR1/GCSH"                  
# [4] "MTHFS/ECHS1/RIDA/PPARA/DBT/ALDH5A1/AUH/DECR1/GCSH"                  
# [5] "SPP1/PLPPR4/AFG3L2/CDH11/FN1/GLI2/NRDC/EFNB1/NTRK2"                 
# [6] "AMPD2/SLC25A1/PDHX/PRPS2/PPARA/ATP6/PARP10/ND2"                     
# [7] "AMPD2/SLC25A1/PDHX/PRPS2/PPARA/ATP6/PARP10/ND2"                     
# [8] "AMPD2/SLC25A1/PDHX/PRPS2/PPARA/ATP6/PARP10/ND2"                     
# [9] "AMPD2/SLC25A1/PDHX/PRPS2/PPARA/ATP6/PARP10/ND2"                     
# [10] "SPP1/ECHS1/PLAAT4/GDPD3/RAB7A/PPARA/AUH/DECR1"  

x <- strsplit("MTHFS/APOBEC3H/ECHS1/RIDA/PPARA/APOBEC3G/DBT/ALDH5A1/AUH/DECR1/GCSH", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1     MTHFS ENSG00000136371
# 2  APOBEC3H ENSG00000100298
# 3     ECHS1 ENSG00000127884
# 4      RIDA ENSG00000132541
# 5     PPARA ENSG00000186951
# 6  APOBEC3G ENSG00000239713
# 7       DBT ENSG00000137992
# 8   ALDH5A1 ENSG00000112294
# 9       AUH ENSG00000148090
# 10    DECR1 ENSG00000104325
# 11     GCSH ENSG00000140905


##------DESeq2_UVM-----
method.name = "DESeq2"
model.type = "Uni"
model.type.2 = "UVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
DESeq2_UVM.cDEGs <- intersect(data.env$DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                              data.env$DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
DESeq2_UVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                    method.name, model.type.2, var.name))
DESeq2_UVM.unique.cDEGs <- DESeq2_UVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
DESeq2_UVM.cDEGs.8809678 <- sub("\\..*", "", DESeq2_UVM.cDEGs) 
DESeq2_UVM.unique.cDEGs.8809678 <- sub("\\..*", "", DESeq2_UVM.unique.cDEGs) 

# GO over-representation analysis
DESeq2_UVM.go.ora <-
  enrichGO(
    gene          = DESeq2_UVM.unique.cDEGs.8809678,
    # universe      = DESeq2_UVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
DESeq2_UVM.go.ora@result |> View()  
# DESeq2_UVM.go.ora.dt <- as.data.table(DESeq2_UVM.go.ora)
# clusterProfiler::dotplot(object = DESeq2_UVM.go.ora, 
#                          color = "pvalue",
#                          showCategory = 10) +
#   theme(axis.text.y = element_text(size = 7))


DESeq2_UVM.go.ora@result[order(DESeq2_UVM.go.ora@result$Count, decreasing = T), ] |> View()
DESeq2_UVM.go.ora@result[order(DESeq2_UVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "epidermis development"                          
# [2] "keratinocyte differentiation"                   
# [3] "epidermal cell differentiation"                 
# [4] "skin development"                               
# [5] "intermediate filament organization"             
# [6] "intermediate filament cytoskeleton organization"
# [7] "intermediate filament-based process"            
# [8] "keratinization"                                 
# [9] "positive regulation of mitotic cell cycle"      
# [10] "connective tissue development"      
DESeq2_UVM.go.ora@result[order(DESeq2_UVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "KRT6B/KRT6A/KRT14/KRT5/COL17A1/KRT77/ZBED2" "KRT6B/KRT6A/KRT14/KRT5/KRT77/ZBED2"        
# [3] "KRT6B/KRT6A/KRT14/KRT5/KRT77/ZBED2"         "KRT6B/KRT6A/KRT14/KRT5/KRT77/ZBED2"        
# [5] "KRT6B/KRT6A/KRT14/KRT5/KRT77"               "KRT6B/KRT6A/KRT14/KRT5/KRT77"              
# [7] "KRT6B/KRT6A/KRT14/KRT5/KRT77"               "KRT6B/KRT6A/KRT5/KRT77"                    
# [9] "CDC25C/UBE2C/CDC7/FOXA1"                    "BARX2/HMGCS2/FOXA1/OSR2" 
x <- strsplit("KRT6B/KRT6A/KRT14/KRT5/COL17A1/KRT77/ZBED2", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1   KRT6B ENSG00000185479
# 2   KRT6A ENSG00000205420
# 3   KRT14 ENSG00000186847
# 4    KRT5 ENSG00000186081
# 5 COL17A1 ENSG00000065618
# 6   KRT77 ENSG00000189182
# 7   ZBED2 ENSG00000177494


##------DESeq2_MVM-----
method.name = "DESeq2"
model.type = "Multi"
model.type.2 = "MVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
DESeq2_MVM.cDEGs <- intersect(data.env$DESeq2.eFDR.GeneName.train[[var.name]][1:threshold], 
                              data.env$DESeq2.eFDR.GeneName.test[[var.name]][1:threshold])
DESeq2_MVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                    method.name, model.type.2, var.name))
DESeq2_MVM.unique.cDEGs <- DESeq2_MVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
DESeq2_MVM.cDEGs.8809678 <- sub("\\..*", "", DESeq2_MVM.cDEGs) 
DESeq2_MVM.unique.cDEGs.8809678 <- sub("\\..*", "", DESeq2_MVM.unique.cDEGs) 

# GO over-representation analysis
DESeq2_MVM.go.ora <-
  enrichGO(
    gene          = DESeq2_MVM.unique.cDEGs.8809678,
    # universe      = DESeq2_MVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
DESeq2_MVM.go.ora@result |> View()  

DESeq2_MVM.go.ora@result[order(DESeq2_MVM.go.ora@result$Count, decreasing = T), ] |> View()
DESeq2_MVM.go.ora@result[order(DESeq2_MVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "response to alcohol"                    "cellular response to alcohol"          
# [3] "cellular response to salt"              "response to salt"                      
# [5] "cAMP biosynthetic process"              "meiotic sister chromatid segregation"  
# [7] "meiotic sister chromatid cohesion"      "meiosis II"                            
# [9] "meiosis II cell cycle process"          "cyclic nucleotide biosynthetic process"
DESeq2_MVM.go.ora@result[order(DESeq2_MVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "ADCY8/CES1/CSF3/ADCY7" "ADCY8/CES1/ADCY7"      "ADCY8/ACHE/ADCY7"     
# [4] "ADCY8/ACHE/ADCY7"      "ADCY8/ADCY7"           "BUB1/HORMAD2"         
# [7] "BUB1/HORMAD2"          "BUB1/HORMAD2"          "BUB1/HORMAD2"         
# [10] "ADCY8/ADCY7"          
x <- strsplit("ADCY8/CES1/CSF3/ADCY7", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1  ADCY8 ENSG00000155897
# 2   CES1 ENSG00000198848
# 3   CES1 ENSG00000262243
# 4   CSF3 ENSG00000108342
# 5  ADCY7 ENSG00000121281


##------edgeR_UVM-----
method.name = "edgeR"
model.type = "Uni"
model.type.2 = "UVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
edgeR_UVM.cDEGs <- intersect(data.env$edgeR.eFDR.GeneName.train[[var.name]][1:threshold], 
                              data.env$edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
edgeR_UVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                    method.name, model.type.2, var.name))
edgeR_UVM.unique.cDEGs <- edgeR_UVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
edgeR_UVM.cDEGs.8809678 <- sub("\\..*", "", edgeR_UVM.cDEGs) 
edgeR_UVM.unique.cDEGs.8809678 <- sub("\\..*", "", edgeR_UVM.unique.cDEGs) 

# GO over-representation analysis
edgeR_UVM.go.ora <-
  enrichGO(
    gene          = edgeR_UVM.unique.cDEGs.8809678,
    # universe      = edgeR_UVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
edgeR_UVM.go.ora@result |> View()  

edgeR_UVM.go.ora@result[order(edgeR_UVM.go.ora@result$Count, decreasing = T), ] |> View()
edgeR_UVM.go.ora@result[order(edgeR_UVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "epidermis development"              "skin development"                  
# [3] "hormone secretion"                  "hormone transport"                 
# [5] "lipid catabolic process"            "lipid transport"                   
# [7] "skin epidermis development"         "regulation of hormone secretion"   
# [9] "organic hydroxy compound transport" "amide transport"     
edgeR_UVM.go.ora@result[order(edgeR_UVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "KRT25/KRT71/SPRR1A/WNT10A/FLG2/FOXQ1/KRTDAP/ALOXE3/LIPN/KLK7/ABCA12/C1orf68/LCE1A/LCE1C"
# [2] "KRT25/KRT71/SPRR1A/WNT10A/FLG2/FOXQ1/ALOXE3/LIPN/ABCA12/LCE1A/LCE1C"                    
# [3] "PTPRN/CHGA/POMC/CELA2A/PLA2G3/ABCA12/FAM3D/F2"                                          
# [4] "PTPRN/CHGA/POMC/CELA2A/PLA2G3/ABCA12/FAM3D/F2"                                          
# [5] "CYP2W1/CEL/PNLIP/CLPS/APOC3/PLA2G1B/LIPH/LIPN"                                          
# [6] "CEL/POMC/PNLIP/APOC3/PLA2G1B/PLA2G3/ABCA12/STAR"                                        
# [7] "KRT25/KRT71/WNT10A/FLG2/FOXQ1/ALOXE3/ABCA12"                                            
# [8] "CHGA/POMC/CELA2A/PLA2G3/ABCA12/FAM3D/F2"                                                
# [9] "CHGA/CEL/POMC/PNLIP/APOC3/ABCA12/STAR"                                                  
# [10] "PTPRN/CHGA/CELA2A/KMO/ABCA12/FAM3D/F2"       
x <- strsplit("KRT25/KRT71/SPRR1A/WNT10A/FLG2/FOXQ1/KRTDAP/ALOXE3/LIPN/KLK7/ABCA12/C1orf68/LCE1A/LCE1C", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1    KRT25 ENSG00000204897
# 2    KRT71 ENSG00000139648
# 3   SPRR1A ENSG00000169474
# 4   WNT10A ENSG00000135925
# 5     FLG2 ENSG00000143520
# 6    FOXQ1 ENSG00000164379
# 7   KRTDAP ENSG00000188508
# 8   ALOXE3 ENSG00000179148
# 9     LIPN ENSG00000204020
# 10    KLK7 ENSG00000169035
# 11  ABCA12 ENSG00000144452
# 12 C1orf68 ENSG00000198854
# 13   LCE1A ENSG00000186844
# 14   LCE1C ENSG00000197084


##------edgeR_MVM-----
method.name = "edgeR"
model.type = "Multi"
model.type.2 = "MVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
edgeR_MVM.cDEGs <- intersect(data.env$edgeR.eFDR.GeneName.train[[var.name]][1:threshold], 
                              data.env$edgeR.eFDR.GeneName.test[[var.name]][1:threshold])
edgeR_MVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                    method.name, model.type.2, var.name))
edgeR_MVM.unique.cDEGs <- edgeR_MVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
edgeR_MVM.cDEGs.8809678 <- sub("\\..*", "", edgeR_MVM.cDEGs) 
edgeR_MVM.unique.cDEGs.8809678 <- sub("\\..*", "", edgeR_MVM.unique.cDEGs) 

# GO over-representation analysis
edgeR_MVM.go.ora <-
  enrichGO(
    gene          = edgeR_MVM.unique.cDEGs.8809678,
    # universe      = edgeR_MVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
edgeR_MVM.go.ora@result |> View()  

edgeR_MVM.go.ora@result[order(edgeR_MVM.go.ora@result$Count, decreasing = T), ] |> View()
edgeR_MVM.go.ora@result[order(edgeR_MVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "maturation of SSU-rRNA"                        
# [2] "ribosomal small subunit biogenesis"            
# [3] "rRNA processing"                               
# [4] "glial cell differentiation"                    
# [5] "ribonucleoprotein complex assembly"            
# [6] "ribonucleoprotein complex subunit organization"
# [7] "rRNA metabolic process"                        
# [8] "positive regulation of growth"                 
# [9] "ribosome biogenesis"                           
# [10] "gliogenesis"       
edgeR_MVM.go.ora@result[order(edgeR_MVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "RPS16/SNU13"  "RPS16/SNU13"  "RPS16/SNU13"  "POU3F2/FGF5" 
# [5] "SNU13/SF3B6"  "SNU13/SF3B6"  "RPS16/SNU13"  "KRT17/POU3F2"
# [9] "RPS16/SNU13"  "POU3F2/FGF5"          
x <- strsplit("RPS16/SNU13", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1  RPS16 ENSG00000105193
# 2  SNU13 ENSG00000100138


##------eBayes_UVM-----
method.name = "eBayes"
model.type = "Uni"
model.type.2 = "UVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
eBayes_UVM.cDEGs <- intersect(names(data.env$eBayes.eFDR.train2[[var.name]][1:threshold]), 
                              names(data.env$eBayes.eFDR.test2[[var.name]][1:threshold]))
eBayes_UVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                   method.name, model.type.2, var.name))
eBayes_UVM.unique.cDEGs <- eBayes_UVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
eBayes_UVM.cDEGs.8809678 <- sub("\\..*", "", eBayes_UVM.cDEGs) 
eBayes_UVM.unique.cDEGs.8809678 <- sub("\\..*", "", eBayes_UVM.unique.cDEGs) 

# GO over-representation analysis
eBayes_UVM.go.ora <-
  enrichGO(
    gene          = eBayes_UVM.unique.cDEGs.8809678,
    # universe      = eBayes_UVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
eBayes_UVM.go.ora@result |> View()  

eBayes_UVM.go.ora@result[order(eBayes_UVM.go.ora@result$Count, decreasing = T), ] |> View()
eBayes_UVM.go.ora@result[order(eBayes_UVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "regulation of fear response"                
# [2] "response to redox state"                    
# [3] "positive regulation of behavior"            
# [4] "behavioral fear response"                   
# [5] "behavioral defense response"                
# [6] "fear response"                              
# [7] "circadian regulation of gene expression"    
# [8] "regulation of behavior"                     
# [9] "multicellular organismal response to stress"
# [10] "positive regulation of DNA repair"    
eBayes_UVM.go.ora@result[order(eBayes_UVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "NPAS2" "NPAS2" "NPAS2" "NPAS2" "NPAS2" "NPAS2" "NPAS2" "NPAS2" "NPAS2"
# [10] "NPAS2"       
x <- strsplit("NPAS2", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1  NPAS2 ENSG00000170485


##------eBayes_MVM-----
method.name = "eBayes"
model.type = "Multi"
model.type.2 = "MVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
eBayes_MVM.cDEGs <- intersect(names(data.env$eBayes.eFDR.train2[[var.name]][1:threshold]), 
                              names(data.env$eBayes.eFDR.test2[[var.name]][1:threshold]))
eBayes_MVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                    method.name, model.type.2, var.name))
eBayes_MVM.unique.cDEGs <- eBayes_MVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
eBayes_MVM.cDEGs.8809678 <- sub("\\..*", "", eBayes_MVM.cDEGs) 
eBayes_MVM.unique.cDEGs.8809678 <- sub("\\..*", "", eBayes_MVM.unique.cDEGs) 

# GO over-representation analysis
eBayes_MVM.go.ora <-
  enrichGO(
    gene          = eBayes_MVM.unique.cDEGs.8809678,
    # universe      = eBayes_MVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
eBayes_MVM.go.ora@result |> View()  

eBayes_MVM.go.ora@result[order(eBayes_MVM.go.ora@result$Count, decreasing = T), ] |> View()
eBayes_MVM.go.ora@result[order(eBayes_MVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "positive regulation of protein kinase activity"                                  
# [2] "positive regulation of kinase activity"                                          
# [3] "positive regulation of cell development"                                         
# [4] "positive regulation of synapse maturation"                                       
# [5] "regulation of T-helper 2 cell differentiation"                                   
# [6] "motor neuron migration"                                                          
# [7] "antigen processing and presentation of exogenous peptide antigen via MHC class I"
# [8] "interneuron migration"                                                           
# [9] "CD40 signaling pathway"                                                          
# [10] "neurotransmitter-gated ion channel clustering"   
eBayes_MVM.go.ora@result[order(eBayes_MVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "CD86/RELN" "CD86/RELN" "CD86/RELN" "RELN"      "CD86"      "RELN"     
# [7] "IFI30"     "RELN"      "CD86"      "RELN" 
x <- strsplit("CD86/RELN", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1   CD86 ENSG00000114013
# 2   RELN ENSG00000189056


##------voom.limma_UVM-----
method.name = "VoomLimma"
method.name.2 = "voom_limma"
model.type = "Uni"
model.type.2 = "UVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
voom.limma_UVM.cDEGs <- intersect(names(data.env$voom.eFDR.train2[[var.name]][1:threshold]), 
                                  names(data.env$voom.eFDR.test2[[var.name]][1:threshold]))
voom.limma_UVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                        method.name.2, model.type.2, var.name))
voom.limma_UVM.unique.cDEGs <- voom.limma_UVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
voom.limma_UVM.cDEGs.8809678 <- sub("\\..*", "", voom.limma_UVM.cDEGs) 
voom.limma_UVM.unique.cDEGs.8809678 <- sub("\\..*", "", voom.limma_UVM.unique.cDEGs) 

# GO over-representation analysis
voom.limma_UVM.go.ora <-
  enrichGO(
    gene          = voom.limma_UVM.unique.cDEGs.8809678,
    # universe      = voom.limma_UVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
voom.limma_UVM.go.ora@result |> View()  

voom.limma_UVM.go.ora@result[order(voom.limma_UVM.go.ora@result$Count, decreasing = T), ] |> View()
voom.limma_UVM.go.ora@result[order(voom.limma_UVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "establishment of protein localization to organelle"             
# [2] "protein insertion into mitochondrial inner membrane"            
# [3] "3'-UTR-mediated mRNA destabilization"                           
# [4] "negative regulation by host of viral transcription"             
# [5] "3'-UTR-mediated mRNA stabilization"                             
# [6] "protein insertion into mitochondrial membrane"                  
# [7] "positive regulation of protein import into nucleus"             
# [8] "establishment of protein localization to mitochondrial membrane"
# [9] "inner mitochondrial membrane organization"                      
# [10] "nuclear membrane organization"     
voom.limma_UVM.go.ora@result[order(voom.limma_UVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "TARDBP/TIMM8A" "TIMM8A"        "TARDBP"        "TARDBP"       
# [5] "TARDBP"        "TIMM8A"        "TARDBP"        "TIMM8A"       
# [9] "TIMM8A"        "TARDBP"    
x <- strsplit("TARDBP/TIMM8A", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1 TARDBP ENSG00000120948
# 2 TIMM8A ENSG00000126953


##------voom.limma_MVM-----
method.name = "VoomLimma"
method.name.2 = "voom_limma"
model.type = "Multi"
model.type.2 = "MVM"
file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/%sModel/Trim/%s%s%s.RData", 
                     model.type, method.name, model.type, seed.i) 
data.env <- local({load(file.name); environment()})
voom.limma_MVM.cDEGs <- intersect(names(data.env$voom.eFDR.train2[[var.name]][1:threshold]), 
                                  names(data.env$voom.eFDR.test2[[var.name]][1:threshold]))
voom.limma_MVM.unique.cDEGs.all.seed <- readRDS(sprintf("../ApplicationResult/UniqueGene/TMM_Top5000/%s_%s.uniq.cDEGs.%s.all.seed.RDS",
                                                        method.name.2, model.type.2, var.name))
voom.limma_MVM.unique.cDEGs <- voom.limma_MVM.unique.cDEGs.all.seed[[as.character(seed.i)]]

# Delete any character after the comma
voom.limma_MVM.cDEGs.8809678 <- sub("\\..*", "", voom.limma_MVM.cDEGs) 
voom.limma_MVM.unique.cDEGs.8809678 <- sub("\\..*", "", voom.limma_MVM.unique.cDEGs) 

# GO over-representation analysis
voom.limma_MVM.go.ora <-
  enrichGO(
    gene          = voom.limma_MVM.unique.cDEGs.8809678,
    # universe      = voom.limma_MVM.cDEGs.8809678,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
voom.limma_MVM.go.ora@result |> View()  

voom.limma_MVM.go.ora@result[order(voom.limma_MVM.go.ora@result$Count, decreasing = T), ] |> View()
voom.limma_MVM.go.ora@result[order(voom.limma_MVM.go.ora@result$Count, decreasing = T), "Description"][1:10]
# [1] "ribonucleoprotein complex assembly"            
# [2] "ribonucleoprotein complex subunit organization"
# [3] "pre-mRNA cleavage required for polyadenylation"
# [4] "mRNA cleavage involved in mRNA processing"     
# [5] "regulation of mRNA polyadenylation"            
# [6] "protein heterotetramerization"                 
# [7] "toll-like receptor 9 signaling pathway"        
# [8] "U2-type prespliceosome assembly"               
# [9] "mRNA cleavage"                                 
# [10] "regulation of mRNA 3'-end processing"   
voom.limma_MVM.go.ora@result[order(voom.limma_MVM.go.ora@result$Count, decreasing = T), "geneID"][1:10]
# [1] "CPSF7/BUD13" "CPSF7/BUD13" "CPSF7"       "CPSF7"       "CPSF7"      
# [6] "CPSF7"       "EPG5"        "BUD13"       "CPSF7"       "CPSF7"  
x <- strsplit("CPSF7/BUD13", split = "/") |> unlist()
suppressMessages({
  ids <- bitr(x, fromType = "SYMBOL", 
              toType = c("ENSEMBL"), 
              OrgDb = "org.Hs.eg.db")
})
ids
# SYMBOL         ENSEMBL
# 1  CPSF7 ENSG00000149532
# 2  BUD13 ENSG00000137656
