# Environment Variables
# Ref: https://cran.r-project.org/web/packages/options/vignettes/envvars.html

file.name <- sprintf("../ApplicationData/derived/RandomSeed/Top5000/UniModel/Trim/DESeq2Uni%s.RData", seed.i)
data.env <- local({load(file.name); environment()})
