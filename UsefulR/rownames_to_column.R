# Convert rownames to the first column and assign the 'La' as the column name
df <- tibble::rownames_to_column(df, "La")
