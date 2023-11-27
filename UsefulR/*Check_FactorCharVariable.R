##------Check the structure and distribution of factor/character variables------
all.fac <- sapply(names(dat.work)[sapply(dat.work, is.factor)], function(x) with(dat.work, table(get(x), useNA = "ifany")), simplify = F)  
View(all.fac)
# Note: variable `Surgery` only have one factor level with non-zero patients.

all.fac.lev <- sapply(names(all.fac), function(x) with(dat.work, levels(get(x))))  
View(all.fac.lev)

all.char <- sapply(names(dat.work)[sapply(dat.work, is.character)], function(x) with(dat.work, table(get(x), useNA = "ifany")))  
View(all.char)
# Note: variable `Surgery2` only have one non-NA group.

##------Set all variables------
vars.all <- c(
  "Site",
  "Race2",
  "Race3",
  "NDI2",
  "Insurance2",
  "BMI",
  "Age.at.Diagnosis",
  "DetectionMode",
  "BRCA.1.2.Status",
  "Laterality.of.DCIS",
  "Grade2",
  "Grade3",
  "Architecture",
  "ER",
  "PR",
  ##------New variable------
  "ER.PR",
  "HER2",
  # "Surgery", # Not use this variable
  ##------New variable------
  # "Surgery2", # Not use this variable
  "SxRT",
  "MarginPositive",
  "ReExcision",
  "Adjuvant.RT",
  "Adjuvant.endocrine.therapy",
  "Laterality.of.Recurrent.DCIS",
  "Any.new.primary.invasive.breast.cancer",
  "Laterality.of.Invasive.BC",
  "Status.at.Last.Follow.Up",
  "Recur",
  "Recur3yr",
  "Recur5yr",
  "T2Recur",
  "T2RecurLFU",
  "T2LFU"
)

vars.cat <-rep(1, length(vars.all))
vars.cat[vars.all %in% c("NDI2", "BMI", "Age.at.Diagnosis","T2Recur", "T2RecurLFU", "T2LFU")] <- 0
