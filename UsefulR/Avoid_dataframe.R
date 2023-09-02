library(microbenchmark)
# Case 1
microbenchmark(mat[, 1], df[, 1])

Unit: microseconds
     expr  min    lq     mean median    uq    max neval cld
 mat[, 1] 1.14 1.375  1.64743   1.53 1.705  8.950   100  a 
  df[, 1] 8.73 9.220 10.03083   9.40 9.705 58.731   100   b


microbenchmark(mat[1, ], df[1, ])

Unit: microseconds
     expr      min       lq       mean   median       uq      max neval cld
 mat[1, ]    4.520    5.130   12.17284   12.085   17.095   30.920   100  a 
  df[1, ] 4451.333 4528.663 4832.24897 4604.429 4723.460 7832.763   100   b

# Case 2
d <- function() {
  data.frame(
    d1 = sample(1:9, 3, replace = TRUE),
    d2 = sample(1:9, 3, replace = TRUE)
  )
}

m <- function() {
  matrix(sample(1:9, 6, replace = TRUE), ncol = 2)
}

microbenchmark(
 data.frame_solution = d(),
 matrix_solution     = m()
)

Unit: microseconds
                expr     min       lq     mean   median      uq       max neval
 data.frame_solution 126.631 133.1305 165.7105 134.5555 141.451  2111.978   100
     matrix_solution   4.920   5.6800 135.3040   7.9950   8.535 12738.269   100
 cld
   a
   a

# Case 3
rolls <- matrix(sample(1:9, 6, replace = TRUE), ncol = 2)
rolls_df <- data.frame(rolls)

app <- function(x) {
    apply(x, 1, sum)
}

microbenchmark(
    df_sol = app(rolls_df),
    mat_sol = app(rolls)
)

Unit: microseconds
    expr   min     lq     mean median     uq     max neval cld
  df_sol 47.27 48.840 55.20131 49.890 51.025 255.971   100   b
 mat_sol 14.52 15.765 16.53882 16.325 17.045  28.531   100  a
