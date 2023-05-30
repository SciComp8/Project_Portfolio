alter.glm.rownames <- function(model) {
  var.names <- names(coef(model))
  for (v in 1:length(var.names)) {
    i <- var.names[v]
    if (i == "(Intercept)") {
      var.names[v] <- "Intercept"
    } else if (i == "Age.Dx") {
      var.names[v] <- "Age at Diagnosis"
    } else if (grepl("*.[a-z][A-Z]", i)) {
      if (grepl("^Race", i)) {
        var.names[v] <- paste0(substring(i, regexpr("*.[a-z][A-Z]", i) + 2), " (vs Non-Asian)")
      } else {
        var.names[v] <- paste0(substring(i, regexpr("*.[a-z][A-Z]", i) + 2), " (vs Private)")
      }
    } else if (i == "MMGSDNo") {
        var.names[v] <- "Mammo Screen-Detected No (vs Mammo Screen-Detected Yes)"
    } else if (i == "MMGOYes") {
        var.names[v] <- "Mammo-Occult Yes (vs Mammo-Occult No)"
    } else if (i == "MRISDYes") {
        var.names[v] <- "MRI Screen-Detected Yes (vs MRI Screen-Detected No)"
    } else if (grepl("*.2[A-Z]", i)) {
        if (grepl("^Hos", i)) {
            var.names[v] <- "NYPQ (vs WCM)"
        } else {
          var.names[v] <- paste0(substring(i, regexpr("*.2[A-Z]", i) + 2), " (vs Heterogenously dense)")
        }
    } else if (i == "ERNegative") {
        var.names[v] <- "ER Negative (vs ER Positive)"
    }
  }
  return(var.names)
}
