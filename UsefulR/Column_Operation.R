# Substract the odd column from the even column
df_fix <- df[, seq(2, ncol(df), by = 2)] - df[, seq(1, ncol(df), by = 2)]
