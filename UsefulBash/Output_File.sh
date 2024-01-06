#!/bin/bash

fitdir="mdl-fits" 
output_dir="$fitdir/output" 

mkdir -p "$fitdir/$fitoutputdir"

i="$JOBINDEX"
filename_base=".mdl-$i"
output_file="$output_dir/$filename_base"
error_file="$output_dir/.err.$filename_base"

Rscript --vanilla mdl.fit.R $i >"$output_file" 2>"$error_file"
# >: redirect standard output; 2>: redirect standard error

# in mdl.fit.R
# args <- commandArgs(trailingOnly=T) 
# i <- as.numeric(args[1]) 
# dir.name <- "model"
# y.i <- readRDS(sprintf("%s/data/y-%s.rds", dir.name, i))
# print(sprintf("Start building the model %s...", i))
# lbd <-  mean(y.i) * (9 ** seq(-10, 3, length.out = 20)) 
# gene.score <- readRDS(sprintf("%s/data/gene-%s.rds", dir.name, i))
# mdl.i <- cv.glmreg(x = gene.score, y = y.i, family = "negbin",
#                    standardize = T,
#                    alpha = 0, theta = 1,
#                    lambda = lbd,
#                    nfolds = 10, n.cores = 10,
#                    maxit = 20000,
#                    thresh = 1e-5)
# print("Complete building the model")
# print("Calculate the correlation between fitted y and observed y:")
# print(cor(mdl.i$fit$fitted.values, y.i, method = "pearson"))
# 
# saveRDS(mdl.i, sprintf("%s/results/mdl.result-%s.rds", dir.name, i))
