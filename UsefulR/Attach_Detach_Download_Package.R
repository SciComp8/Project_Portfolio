bio_pkgs <- c("GenomicRanges", "Rsamtools", "rtracklayer", "GenomicAlignments", "ShortRead")
gen_pkgs <- c("ggplot2")
for (pkg in c(bio_pkgs, gen_pkgs)) {
  if (!require(pkg, character.only = T, quietly = T) & (pkg %in% bio_pkgs)) {
    BiocManager::install(pkg)
    library(pkg, character.only = T)
  } else if (!require(pkg, character.only = T, quietly = T) & (pkg %in% gen_pkgs)) {
    install.packages(pkg)
    library(pkg, character.only = T)
  } else {
    library(pkg, character.only = T)
  }
}

a = "GenomicRanges"
library(a)
# Error in library(a) : there is no package called ‘a’
library(a, character.only = T)

detach('package:clusterProfiler', unload = T)
