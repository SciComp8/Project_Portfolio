# slice_max()
# Given a dataframe with 3 columns - seed, class, frequency where each seed has multiple classes and each class has its own frequency
# Extract the top 10 most frequent class per seed
BMI.class.freq.top10.per.seed <- BMI.class.freq |> 
  group_by(Seed) |> 
  slice_max(order_by = Frequency, n = 10, with_ties = F) |>
  ungroup()

# slice_tail()
# Select the last row of the data frame
t_res <- lm(rank ~ - 1 + condition, data = df_med)
t_freq <- t_res |> 
   tidy() |>
   mutate(condition = str_remove(term, "condition"), .before = 1) |> 
   select(-term)
t_freq |> slice_tail()
