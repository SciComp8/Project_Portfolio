n <- 10000
# Slow code
growing <- function(n) {
    x <- NULL
    for(i in 1:n)
        x <- c(x, rnorm(1))
    x
}
system.time(res_grow <- growing(n = 10000))


# Fast code
pre_allocate <- function(n) {
    x <- numeric(n) # Pre-allocate
    for(i in 1:n) 
        x[i] <- rnorm(1)
    x
}
system.time(res_allocate <- pre_allocate(n = 10000))
