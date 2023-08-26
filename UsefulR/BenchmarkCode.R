colon <- function(n) 1:n
seq_default <- function(n) seq(1, n)
seq_by <- function(n) seq(1, n, by = 1)
seq_al <- function(n) seq_along(rep(LETTERS, length.out = n))

system.time(colon(1e6))
system.time(seq_default(1e6))
system.time(seq_by(1e6))
system.time(seq_al(1e6))

library(microbenchmark)
n <- 1e6
microbenchmark(colon(n),
               seq_default(n),
               seq_by(n),
               seq_al(n),
               times = 10)
