library(parallel)

# Old version
dat.pheno.train <- readRDS("../ApplicationData/derived/RandomSeed/dat.pheno.train120.RDS")
dat.pheno.test <- readRDS("../ApplicationData/derived/RandomSeed/dat.pheno.test120.RDS")
vars <- c("BMI", "AGE", "SEX", "MHABNWBC", "MHARTHTS", "MHCVD", "MHASTHMA", "MHBCTINF", "MHBLDDND", "MHCOCAINE5", "MHCOPD", "MHDRNKSTS", "SMTSPAX") 

design.train <- mclapply(1:length(vars),
                         function(i) model.matrix(~dat.pheno.train[vars][[i]]),
                         mc.cores = 10)

design.test <- mclapply(1:length(vars),
                        function(i) model.matrix(~dat.pheno.test[vars][[i]]),
                        mc.cores = 10)

# New version
design_mat_i <- function(var_focus, data_set) {
  model.matrix(~., data_set[var_focus])
}

var_focus <- c("BMI", "AGE", "SEX", "MHABNWBC", "MHARTHTS", "MHCVD", "MHASTHMA", "MHBCTINF", "MHBLDDND", "MHCOCAINE5", "MHCOPD", "MHDRNKSTS", "SMTSPAX") 
data_all <- list(dat.pheno.train, dat.pheno.test)
design_train <- mcmapply(design_mat_i, var_focus, data_all[1], SIMPLIFY = FALSE, mc.cores = 10)
design_test <- mcmapply(design_mat_i, var_focus, data_all[2], SIMPLIFY = FALSE, mc.cores = 10)

