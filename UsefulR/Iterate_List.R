data_list <- list(
  data.frame(id = 1:3, value = c(10, 20, 30)),
  data.frame(id = 4:6, value = c(40, 50, 60)),
  data.frame(id = 7:9, value = c(70, 80, 90))
)

result_list <- list()

for (df in data_list) {
  df$value <- df$value * 2
  result_list[[length(result_list) + 1]] <- df
  # 0+1; 1+1; 2+1 
}

for (i in seq_along(result_list)) {
  cat("Modified data frame", i, ":\n")
  print(result_list[[i]])
}
