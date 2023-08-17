# Among cDEGs detected within each SMA approach, estimate the percentage of these cDEGs identified by multi-model BMAseq - all seeds
p.fun <- function(SMA.method = NULL, data.set = NULL) {
  d <- data.set |> 
    filter(.data[[SMA.method]] == 1) |> 
    summarise(n = n()) |> pull()
  n <- data.set |> 
    filter(.data[[SMA.method]] == 1, BMAseq == 1) |> 
    summarise(n = n()) |> pull()
  p <- n/d
  return(p)
}

p.list <- vector(mode = "list", length = 10)
seed.vec <- c(8809678, 98907, 233, 556, 7890, 120, 2390, 778, 666, 99999)
names(p.list) <- seed.vec 
for (z in 1:10) {
  class.freq <- readRDS(file = sprintf("../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/%s_5000_%s.RDS", var.name, seed.vec[z]))
  class.freq <- class.freq[, -ncol(class.freq)]
  class.freq <- class.freq |> mutate(across(.cols = everything(), as.numeric))
  class.freq[, ncol(class.freq) + 1] <- rowSums(class.freq)
  colnames(class.freq)[ncol(class.freq)] <- "row.sum"
  p.vec <- vector(mode = "numeric", length = 8)
  for (i in seq_along(SMA.method.vec)) {
    p <- p.fun(SMA.method = SMA.method.vec[i], data.set = class.freq)
    p.vec[i] <- p
  }
  names(p.vec) <- SMA.method.vec
  p.list[[z]] <- p.vec
}

p.list.to.df <- do.call(rbind, p.list)
p.list.to.df.avg <- p.list.to.df |> 
  as.data.frame() |> 
  summarise(across(.cols = everything(), mean, .names = "mean_{.col}"))
