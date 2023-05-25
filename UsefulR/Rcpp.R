sumR <- function(x) {
  total <- 0
  for (i in seq_along(x)) {
    total <- total + x[i]
  }
  total
}

# C++ equivalent
cppFunction('double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; i++) {
    total += x[i];
  }
  return total;
}')

sum(c(2, 3, 6))
sumR(c(2, 3, 6))
sumC(c(2, 3, 6))

if (!require("microbenchmark", quietly = TRUE)) { install.packages("microbenchmark"); library(microbenchmark) } else { library(microbenchmark) }
x <- runif(1e3)
microbenchmark(
  sum(x),
  sumC(x),
  sumR(x)
)
