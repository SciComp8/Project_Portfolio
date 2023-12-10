library(magrittr)
df <- read.csv("prevend.csv")

lm_no_statin <- df |> filter(Statin == 0) %>% lm(RFFT ~ Age, data = .) # The second pipe operator cannot be |>
lm_statin <- df |> filter(Statin == 1) %>% lm(RFFT ~ Age, data = .)

# Question: are 2 slopes significantly different from each other between group 1 (Statin == 0) and group 2 (Statin == 1)?

# Solution 1: do 2 95% CIs NOT overlap from each other?
confint(lm_no_statin, level = 0.95)
confint(lm_statin, level = 0.95)

# Solution 2: is the interaction term Age*Statin significantly different from 0?
lm_full <- lm(RFFT ~ Age * factor(Statin), data = df)
summary(lm_full)
