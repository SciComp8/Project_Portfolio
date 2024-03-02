# https://ycphs.github.io/openxlsx/articles/Introduction.html
# https://ycphs.github.io/openxlsx/articles/Formatting.html

wb <- createWorkbook()
addWorksheet(wb, "data_A1")
addWorksheet(wb, "data_A2")
addWorksheet(wb, "data_B1")
addWorksheet(wb, "data_B2")
writeData(wb, sheet = 1, data_A1,  startRow = 1, startCol = 1)
writeData(wb, sheet = 1, data_A2, startRow = 1, startCol = 17)
writeData(wb, sheet = 2, data_B1,  startRow = 1, startCol = 1)
writeData(wb, sheet = 2, data_B2, startRow = 1, startCol = 17)
saveWorkbook(wb, "../result/outcome.xlsx", overwrite = TRUE)

# Iterate through the sheet name and data objects contained in a list
dataset_list <- list(
  list(name = "data_A1", data = data_A1),
  list(name = "data_A2", data = data_A2),
  list(name = "data_B1", data = data_B1),
  list(name = "data_B2", data = data_B2)
)
wb <- createWorkbook()
for (i in seq_along(dataset_list)) {
  sheet_name <- dataset_list[[i]]$name
  data <- dataset_list[[i]]$data
  addWorksheet(wb, sheet_name)
  if (i % 2 != 0 ) { 
    writeData(wb, sheet = i, data, startRow = 1, startCol = 1)
  } else {
    writeData(wb, sheet = i/2, data, startRow = 1, startCol = 17)
  }
} 

saveWorkbook(wb, "../result/outcome.xlsx", overwrite = TRUE)
