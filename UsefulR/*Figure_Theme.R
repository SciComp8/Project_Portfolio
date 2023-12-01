# Reference: https://cran.r-project.org/web/packages/ggplot2/index.html

theme_BMA <- function(
    base_size = 12,
    base_family = "Arial"
){
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = "black", vjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6), # Customize the box border in the strip
      # strip.background = element_blank() # This is gonna remove the strip background
      # strip.background = element_rect(fill = "transparent")
      strip.placement = "inside"
    )
}

# Pretty contrast color
# steelblue vs indianred
# #99CCFF (lightblue) vs #FFCCCC (lightpink)
# green vs grey
# lightblue vs lightgreen
# lightblue vs darkblue
# forestgreen vs maroon

# Pretty single color
# olivedrab4
# navyblue
