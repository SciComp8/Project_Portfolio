# Estimate the mean of target_var in each combinatin of group_1, group_2, group_3
# If target_var is binary (0, 1), such mean will be the proportion of being 1 in each combinatin of group_1, group_2, group_3
with(data, tapply(X=target_var, INDEX=list(group_1, group_2, group_3), FUN=mean)) |>
  ftable() |>
  round(2)
