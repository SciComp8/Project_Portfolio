# Factorize a character variable and order its levels based on the count of 1's in each element
count_ones <- function(string) {
  sum(as.integer(strsplit(string, "")[[1]]))
}

uniq_class <- unique(BMI.class.freq.top10.per.seed$Class)

BMI.class.freq.top10.per.seed$Class <- factor(BMI.class.freq.top10.per.seed$Class, levels = uniq_class[order(sapply(uniq_class, count_ones), decreasing = T)])
