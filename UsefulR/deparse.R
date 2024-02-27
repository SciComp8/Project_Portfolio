# Convert unevaluated expressions into character strings
deparse(args(DESeq2::DESeq), width.cutoff = 500)
[1] "function (object, test = c(\"Wald\", \"LRT\"), fitType = c(\"parametric\", \"local\", \"mean\", \"glmGamPoi\"), sfType = c(\"ratio\", \"poscounts\", \"iterate\"), betaPrior, full = design(object), reduced, quiet = FALSE, minReplicatesForReplace = 7, modelMatrixType, useT = FALSE, minmu = if (fitType == \"glmGamPoi\") 1e-06 else 0.5, parallel = FALSE, BPPARAM = bpparam()) "
[2] "NULL"         
