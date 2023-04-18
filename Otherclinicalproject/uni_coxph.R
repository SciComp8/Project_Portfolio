"%_%" <- function(m, n) paste0(m, "_", n)
"%0%" <- function(m, n) paste0(m, n)
uni_coxph <- function(surv.time, surv.event) {
  out <- mclapply(1:length(vars.ana),
                  function(i) {
                    vars.missing <- ifelse(is.na(dat.work[, vars.ana[i]]), 
                                         "missing", "non.missing") %>%
                      factor(levels = c("non.missing", "missing"))
                    res.i <- try(coxph(formula("Surv(" %0% surv.time %0% "," %0% surv.event %0% ") ~ " %0% vars.ana[i]), data = dat.work), silent = T) # Store errors due to insufficient levels of a factor variable
                    if(any(vars.missing == "missing")) {
                      res.i.miss <- coxph(formula("Surv(" %0% surv.time %0% "," %0% surv.event %0% ") ~ vars.missing"), 
                                          data = dat.work)
                      return(list(res.i, res.i.miss))
                    } else {
                      return(list(res.i))
                    }
                  },
                mc.cores = 10)
  names(out) <- vars.ana
  return(out)
}
