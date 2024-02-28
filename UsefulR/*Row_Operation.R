# Move the first row to the last row
df_fix <- rbind(df[-1, ], df[1, ]) 

# Summation per row
library(microbenchmark)
rolls <- matrix(sample(1:9, 6, replace = TRUE), ncol = 2)

app <- function(x) {
    apply(x, 1, sum)
}

r_sum <- function(x) {
    rowSums(x)
}

microbenchmark(
    app_sol = app(rolls),
    r_sum_sol = r_sum(rolls)
)

Unit: microseconds
      expr   min      lq     mean median     uq      max neval cld
   app_sol 14.30 14.9850 27.93317 15.695 16.615 1069.004   100   a
 r_sum_sol  3.17  3.4245 15.91859  3.655  4.130 1172.435   100   a

# Think of DNA methylation analysis scenarior: remove any CpG probes that have any samples showing detection p value >= 0.01
probe.keep <- rowSums(det.p.filter < 0.01) == ncol(norm.set) 

