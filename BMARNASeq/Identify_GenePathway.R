library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
keytypes(org.Hs.eg.db)

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
