# Assume X ~ Pois(6). Find P(X = 2)
dpois(2,6)

# Assume X ~ Pois(6). Find P(X <= 2)
ppois(2,6)

# Assume X ~ Pois(6). Find P(X > 2)
1 - ppois(2,6)

# Assume X ~ Gamma(6, 1/2). Find P(0.5 < X < 1.5)
pgamma(1.5, 6, 1/2) - pgamma(0.5, 6, 1/2)

# Assume X ~ Norm(6, 9). Find y that P(Y < y) = 0.90
qnorm(0.90, 6, 9)

# Assume X ~ Norm(6, 9). Find y that P(-y < Y < y) = 0.90
qnorm(0.95, 6, 9)
