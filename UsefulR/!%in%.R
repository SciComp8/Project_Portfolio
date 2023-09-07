probe.keep <- !(featureNames(norm.set) %in% anndata.450k$Name[anndata.450k$chr %in% c("chrX","chrY")])
