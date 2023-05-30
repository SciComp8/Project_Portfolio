"%0%" <- function(m, n) paste0(m, n)

alter.glm.rownames <- function(model) {
  var.names <- names(coef(model))
  for (v in 1:length(var.names)) {
    i <- var.names[v]
    if (i == "(Intercept)") {
      var.names[v] <- "Intercept"
    } else if (i == "Age.Dx") {
      var.names[v] <- "Age at Diagnosis"
    } else if (grepl("*.[a-z][A-Z]", i)) {
      if (grepl("^Race.A", i)) {
        var.names[v] <- substring(i, regexpr("*.[a-z][A-Z]", i) + 2) %0% " (vs Non-Asian)"
      } else if (grepl("^RaceE", i)) {
        var.names[v] <- substring(i, regexpr("*.[0-9][A-Z]", i) + 2) %0% " (vs Non-Hispanic White)"
      } else if (grepl("^Ins", i)) {
        var.names[v] <- substring(i, regexpr("*.[a-z][A-Z]", i) + 2) %0% " (vs Private)"
      } else if (grepl("^Trip", i)) {
        var.names[v] <- "Triple Negative Yes (vs Triple Negative No)"
      }
    } else if (i == "MMGSDNo") {
      var.names[v] <- "Mammo Screen-Detected No (vs Mammo Screen-Detected Yes)"
    } else if (i == "MMGOYes") {
      var.names[v] <- "Mammo-Occult Yes (vs Mammo-Occult No)"
    } else if (i == "MRISDYes") {
      var.names[v] <- "MRI Screen-Detected Yes (vs MRI Screen-Detected No)"
    } else if (grepl("*.[0-9][A-Z]", i)) {
        if (grepl("^Hos", i)) {
          var.names[v] <- "NYPQ (vs WCM)"
        } else if (grepl("^MMGD", i)) {
          var.names[v] <- substring(i, regexpr("*.[0-9][A-Z]", i) + 2) %0% " (vs Heterogenously dense)"
        } else {
          var.names[v] <- "HER2 Positive (vs HER2 Negative)" }
    } else if (i == "ERNegative") {
      var.names[v] <- "ER Negative (vs ER Positive)"
    } else if (i == "PRNegative") {
      var.names[v] <- "PR Negative (vs PR Positive)"
    }
  }
  return(var.names)
}


# Case study
covariates <- predictors[-which(predictors %in% c("RaceEthnicity.4", "BMI", "PR", "HER2", "Subtype", "TripleNegative"))] 
data.work <- na.omit(data.full[c("StudyID", response, predictors)])

## Run multivariate logistic models using Early.stage, Age.Dx, Race.Asian, Hospital.2, Insure.type, MMGSD, MMGO, MRISD, MMGD.2, ER
glm.fit <- glm("%0%"(response%0%"~", "%+%"(covariates)), data = data.work, family = binomial)

## Output odds ratio and its 95% CI
result.multi <- data.frame(OR = exp(coef(glm.fit)), 
                           lower.limit = exp(confint(glm.fit))[, 1], 
                           upper.limit = exp(confint(glm.fit))[, 2], 
                           p.value = summary(glm.fit)$coef[, 4]) 

## Format the output
rownames(result.multi) <- NULL
result.multi$object <- alter.glm.rownames(glm.fit)

result.fit <- result.multi %>%
  dplyr::select(object, everything()) %>%
  mutate(across(.cols = c(OR, lower.limit, upper.limit), .fns = round, 2)) %>% 
  mutate(p.value = round(p.value, 3)) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "< 0.001", p.value)) %>%
  filter(object != "Intercept")

## Present the final logistic regression modeling results
result.fit |> knitr::kable(col.names = c("Covariate", "Adjusted Odds Ratio", "95% CI Lower Limit", "95% CI Upper Limit", "P-value"), # Add column names
                           align = "lcccc")
