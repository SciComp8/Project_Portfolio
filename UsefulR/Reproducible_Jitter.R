library(ggplot2)
jitter_set <- position_jitter(height = 0.1, seed = 689) # A random seed to make the jitter reproducible
ggplot(
  data = df_med,
  mapping = aes(x = rank, y = condition)) +     
  geom_point(position = jitter_set, alpha = 0.5, color = "#7c20c1") + 
  scale_x_continuous(breaks = min(df_med$rank):max(df_med$rank))
