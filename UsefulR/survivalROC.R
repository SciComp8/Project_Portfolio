# More info: http://faculty.washington.edu/heagerty/Software/SurvROC/SurvivalROC/survivalROCdiscuss.pdf
library(survivalROC)
library(tidyverse)
dat0 <- read.csv("tumor_mutation_burden.csv")

dat0$OS_YEARS <- dat0$OS_MONTHS / 12

dat0 <- dat0 |>
  mutate(OS_STATUS2 = ifelse(dat0$OS_STATUS == "1:DECEASED", 1, 0)) |>
  select(OS_MONTHS, OS_YEARS, OS_STATUS2, CANCER_TYPE, TMB_NONSYNONYMOUS)
  
dat.NSCLC <- dat0 |>
  filter(CANCER_TYPE == "Non-Small Cell Lung Cancer")
  
sROC <- survivalROC(Stime = dat.NSCLC$OS_YEARS, status = dat.NSCLC$OS_STATUS2, marker = dat.NSCLC$TMB_NONSYNONYMOUS, predict.time = 5, method = "NNE")
  
plot(sROC$FP, sROC$TP, type = "l", xlim = c(0,1), ylim = c(0,1), col = "blue", 
     xlab = "False positive rate", ylab = "True positive rate",
     lwd = 2, cex.main = 1.3, cex.lab = 1.5, cex.axis = 1.2, font = 1.2)
abline(0,1)
auc.text = paste0("5 year (AUC = ", sprintf("%.3f", sROC$AUC), " )")
legend("bottomright", c(auc.text), lwd = 2, bty = "n", col = c("blue"))
