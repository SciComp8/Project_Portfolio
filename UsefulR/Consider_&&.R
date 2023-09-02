library(microbenchmark)
is_double <- c(F, T, T)

move <- function(is_double) {
    if (is_double[1] & is_double[2] & is_double[3]) {
        current <- 9
    }
}

improved_move <- function(is_double) {
    if (is_double[1] && is_double[2] && is_double[3]) {
        current <- 9
    }
}

microbenchmark(move(is_double), improved_move(is_double), times = 1e5)

Unit: nanoseconds
                     expr min  lq     mean median  uq     max neval cld
          move(is_double) 469 500 697.7305    510 530 7780649 1e+05   b
 improved_move(is_double) 250 280 362.4099    290 300 2749240 1e+05  a
