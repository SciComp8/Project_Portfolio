probe.keep.1 <- !(featureNames(norm.set) %in% anndata.450k$Name[anndata.450k$chr %in% c("chrX","chrY")])

probe.keep.2 <- !(featureNames(norm.set.filter) %in% xReactiveProbes$norm.set.filter)
