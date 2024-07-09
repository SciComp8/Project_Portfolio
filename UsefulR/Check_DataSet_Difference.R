# Check differences in 2 R dataframe objects
any(df_1 != df_2)

# Check differences in 2 spreadsheets
library(daff)
Microcalcs.v5 <- readxl::read_xlsx(
  "../data/raw/Coded Microcalcs v5 offline.xlsx",
  sheet = "Coded",
  range = c("A1:AC122"),
  col_types = c(
    rep("numeric", 14),
    rep("date", 3),
    rep("numeric", 6),
    "text", # X21: Dorexcision: 3/23/2022; 4/13/22; 5/4/22
    rep("numeric", 2), 
    "date",
    "numeric",
    "text"
  ),
  col_names = TRUE,
  na = c("Unknown", "unknown", "N/A", "n/a")
) 

Microcalcs.v6.old <- readxl::read_xlsx(
  "../data/raw/Coded Microcalcs v6 offline_Old.xlsx",
  sheet = "Coded",
  range = c("A1:AC122"),
  col_types = c(
    rep("numeric", 14),
    rep("date", 3),
    rep("numeric", 6),
    "text", # X21: Dorexcision: 3/23/2022; 4/13/22; 5/4/22
    rep("numeric", 2), 
    "date",
    "numeric",
    "text"
  ),
  col_names = TRUE,
  na = c("Unknown", "unknown", "N/A", "n/a")
) 

Microcalcs.v6.new <- readxl::read_xlsx(
  "../data/raw/Coded Microcalcs v6 offline.xlsx",
  sheet = "Coded",
  range = c("A1:AC122"),
  col_types = c(
    rep("numeric", 14),
    rep("date", 3),
    rep("numeric", 6),
    "text", # X21: Dorexcision: 3/23/2022; 4/13/22; 5/4/22
    rep("numeric", 2), 
    "date",
    "numeric",
    "text"
  ),
  col_names = TRUE,
  na = c("Unknown", "unknown", "N/A", "n/a")
) 

patch <- diff_data(Microcalcs.v6.old, Microcalcs.v6.new)
# write_diff(patch, "patch.csv")
render_diff(patch)
