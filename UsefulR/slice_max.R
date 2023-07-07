# Give a dataframe with 3 columns - seed, class, frequency where each seed has multiple classes and each class has its own frequency
# Extract the top 10 most frequent class per seed
BMI.class.freq.top10.per.seed <- BMI.class.freq |> 
  group_by(Seed) |> 
  slice_max(order_by = Frequency, n = 10, with_ties = F) |>
  ungroup()
