# Save and load image
image.date <- format(Sys.Date(), "%Y%b%d")
save.image(file = paste0("../ApplicationData/derived/", image.date, "_modelspace_image_ALiu.RData"))

load("../ApplicationData/derived/2023Mar11_modelspace_image_ALiu.RData")


# Attach libraries and functions
suppressPackageStartupMessages(easypackages::libraries("BMAseq", "limma", "qvalue", "parallel", "ggVennDiagram", "gridExtra", "tidyverse", "edgeR", "DESeq2", "microbenchmark", "gtsummary", "fst"))
source("./2022Dec8Version_ALiu/R/BMAseq.multi.postprob.norm.R")
source("./2022Dec8Version_ALiu/R/Bayesfactor.R")


# Load the preprocessed data
seed.vec.trial <- c(8809678)
for (seed.i in seed.vec.trial) {
  dat.expr.train <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.train%s.RDS", seed.i))
  dat.expr.test <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.expr.test%s.RDS", seed.i))
  dat.pheno.train <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.train%s.RDS", seed.i))
  dat.pheno.test <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/dat.pheno.test%s.RDS", seed.i))
}


## Perform the multivariate Analysis with Interaction
```{r}
var.pool <- c("BMI", "AGE", "SEX", "MHABNWBC") 

var.pool0 <- c("BMI", "SEX") 
MS.default.train0 <- Modelspace(dat.pheno = dat.pheno.train, 
                                var.pool = var.pool0, 
                                max.nvar = 2, 
                                interaction = "BMI&SEX")
myModelSpace <- c(MS.default.train0$model.space,
                  paste0(myModelSpace0, "+AGE"),
                  paste0(myModelSpace0, "+MHABNWBC"),
                  paste0(myModelSpace0, "+AGE+MHABNWBC"))

MS.default.train <- Modelspace(dat.pheno = dat.pheno.train, 
                               var.pool = var.pool, # var.pool should not contain "BMI&SEX"
                               max.nvar = 4, # Maximum number of main-effect variables 
                               interaction = "BMI&SEX")
# MS.default.train$dat.pheno.new

output.multi.int1.train <- BMAseq.multi.postprob.MSout(
  dat.expr.counts = dat.expr.train, 
  dat.pheno = MS.default.train$dat.pheno.new,
  model.space = myModelSpace, 
  cut.BF = 1)

output.multi.int2.train <- BMAseq.multi.DEG(
  dat.pheno = output.multi.int1.train$dat.pheno.new,
  model.space = output.multi.int1.train$model.space, 
  post.modelprob = output.multi.int1.train$post.modelprob, 
  var.pool = var.pool,
  interact = T, 
  cut.FDR = 0.05)
```
