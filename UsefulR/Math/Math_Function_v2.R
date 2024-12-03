# Calculate the sensitivity/specificity
# sensitivity = no. of true positive / (no. of true positive + no. of false negative)
## e.g., sensitivity is the probability of a positive treatment response given that the patient has been treatment successfully.

# specificity = no. of true negative / (no. of true negative + no. of false positive)
## e.g., specificity is the probability of a negative treatment response given that the patient has not been treatment successfully.

# positive predictive value = (no. of true positives) / (no. of true positives + no. of false positives) 
## e.g., positive predictive value is the probability of being treatmented successfully given a positive treatment response.

# negative predictive value = (no. of true negatives) / (no. of true negatives + no. of false negatives) 
## e.g., negative predictive value is the probability of not being treatmented successfully given a negative treatment response.


# Global Sensitivity Analysis of Model Outputs
https://cran.r-project.org/web/packages/sensitivity/index.html

# The replicate function is used to perform bootstrapping B times.
# In each iteration, it samples n observations with replacement from the 'var' variable and calculates the proportion of "Yes" values in the sample.
# The results are stored in the bootstrap_prop vector, which contains the proportions for each bootstrap sample.

set.seed(008)
B <- 200
n <- nrow(df)
bootstrap_prop <- replicate(B, {
  x <- sample(df$var, size = n, replace = TRUE)
  sum(x == "Yes") / n
})

library(boot)
data$age_group <- ifelse(data$AGE >= 20 & data$AGE <= 30, 1, 0)
calc_proportion <- function(data, indices) {
  data_resample <- data[indices, ]
  return(sum(data_resample$age_group) / nrow(data_resample))
}
results <- boot(data = data, statistic = calc_prop, 500)
