library(microbenchmark)
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
