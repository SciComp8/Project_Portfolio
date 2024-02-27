my_df <- data.frame(A = 3:5, B = c("x", "y", "z"))

df_name <- deparse(substitute(my_df))

print(df_name)
