# Users are able to input a data.table object to execute the following function

cor_plot_wrapper <-
  function(cor_method = "spearman",
           var_y = "BMI",
           y_label = "Body mass index",
           var_focus = NULL,
           data_set = NULL) {
    yname <- var_y
    ylab <- y_label
    y <- as.vector(data_set[, get(yname)]) # Modify
    vars.ana <- var_focus
    vars.ana.class <-
      unlist(lapply(data_set[, mget(vars.ana)], function(x)
        paste(class(x), collapse = "."))) # Modify
    vars.ana.ch <- vars.ana[grep("factor", vars.ana.class)]
    
    vars.ana.cat <- rep(0, length(vars.ana))
    vars.ana.cat[vars.ana %in% vars.ana.ch] <- 1
    
    out <- vector(length = length(vars.ana), mode = "list")
    names(out) <- vars.ana
    
    for (i in 1:length(vars.ana)) {
      cat(paste0("\n#### ", yname, " ~ ", vars.ana[i], "\n"))
      x <- data_set[, get(vars.ana[i])] # Modify
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
          y.plab = ylab,
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
          y.plab = ylab,
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
  }
