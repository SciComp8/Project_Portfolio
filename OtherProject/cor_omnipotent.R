easypackages::libraries("tidyverse", "readxl", "gtsummary", "flextable", "stringi", "fst", "broom", "labelled", "data.table", "BTKR", "corrr", "multcomp", "survival", "GGally") |> suppressPackageStartupMessages() |> suppressPackageStartupMessages()

cor_omnipotent <- function(var_y = "BMI", 
                           y_label = "Body mass index",
                           data_set = data_set,
                           var_focus = c("BMI",
                                         "Height",
                                         "Weight",
                                         "Tumor.size.cm",
                                         "Oncotype.score",
                                         "Age.at.diagnosis"),
                           cor_method = "spearman",
                           cor_format = c("pairwise_cor_plot", "cor_table", "ggpairs")) {
  cor_format <-
    match.arg(cor_format, c("pairwise_cor_plot", "cor_table", "ggpairs"))
  if (cor_format == "pairwise_cor_plot") {
    y <- data_set[, var_y]
    vars.ana <- var_focus[-which(var_y %in% var_focus)]
    vars.ana.class <-
      unlist(lapply(data_set[, vars.ana], function(x)
        paste(class(x), collapse = ".")))
    vars.ana.ch <- vars.ana[grep("factor", vars.ana.class)]
    
    vars.ana.cat <- rep(0, length(vars.ana))
    vars.ana.cat[vars.ana %in% vars.ana.ch] <- 1
    
    out <- vector(length = length(vars.ana), mode = "list")
    names(out) <- vars.ana
    
    for (i in 1:length(vars.ana)) {
      cat(paste0("\n#### ", var_y, " ~ ", vars.ana[i], "\n"))
      x <- data_set[, vars.ana[i]]
      if (vars.ana.cat[i] == 1) {
        xlevs <- levels(factor(x))
        xlevs <- xlevs[xlevs %in% unique(x)]
        xnlevs <- length(xlevs)
        
        if (xnlevs == 2) {
          out[[i]] <- wilcox.test(y ~ x, data = data_set)
          p <- out[[i]]$p.value
          cat("\n~~~~~\n")
          print(out[[i]])
          cat("\n~~~~~\n")
          compare_method <- "wilcox"
        } else {
          out[[i]] <- kruskal.test(y ~ x, data = data_set)
          p <- out[[i]]$p.value
          cat("\n~~~~~\n")
          print(out[[i]])
          cat("\n~~~~~\n")
          compare_method <- "kruskal"
        }
        
        p.txt <- fpval.txt(p)
        p.txt <- ifelse(p.txt == "<0.001", "P<0.001", p.txt)
        
        fsmry.graph(
          y = y,
          x = factor(as.vector(x), levels = xlevs),
          type = "bxp",
          y.plab = y_label,
          x.plab = label_attribute(data_set[[vars.ana[i]]]),
          xnames = xlevs,
          mar = c(3.5, 3.5, 0.5, 0.5),
          mgp = c(2.0, 0.8, 0),
          cex.lab = 1,
          cex.axis = 1
        )
        legend(
          "topleft",
          legend = paste0(compare_method, ":", p.txt),
          bty = "n",
          cex = 1
        )
        cat("\n")
      } else {
        out[[i]] <- cor.test(y, x, method = cor_method)
        cat("\n~~~~~\n")
        print(out[[i]])
        cat("\n~~~~~\n")
        
        p <- out[[i]]$p.value
        p.txt <- fpval.txt(p)
        p.txt <- ifelse(p.txt == "<0.001", "P<0.001", p.txt)
        rho <- round(out[[i]]$estimate, 2)
        fsmry.graph(
          y = y,
          x = x,
          type = "scatter",
          y.plab = y_label,
          x.plab = label_attribute(data_set[[vars.ana[i]]]),
          mar = c(3.5, 3.5, 0.5, 0.5),
          mgp = c(2.0, 0.8, 0),
          cex.lab = 1,
          cex.axis = 1
        )
        xytline <- line(x = x, y = y)
        abline(coef(xytline))
        legend(
          "topleft",
          paste0(cor_method, " r=", rho, ", ", p.txt),
          cex = 1.0,
          bty = "n"
        )
        cat("\n")
      }
    }
    
  } else if (cor_format == "cor_table") {
    x <- data_set |>
      dplyr::select(all_of(var_focus)) |>
      correlate(
        use = "pairwise.complete.obs",
        method = cor_method,
        diagonal = 1,
        quiet = T
      ) |>
      shave()
    fashion(x) |> knitr::kable()
  } else {
    var_idx <- match(var_focus, names(data_set))
    g <- ggpairs(
      data = data_set,
      columns = var_idx,
      xlab = "Variables",
      ylab = "Variables",
      title = "Relationships between variables",
      upper = list(continuous = wrap("cor", method = cor_method))
    ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(g)
  }
}
