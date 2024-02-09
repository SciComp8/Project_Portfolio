numGroups <- 10
numRepeats <- 100
totalObservations <- numGroups * numRepeats

groupFactor <- factor(round(numGroups * runif(totalObservations)))
data <- rnorm(totalObservations) + sqrt(as.numeric(groupFactor))
groupedData <- split(data, groupFactor)
