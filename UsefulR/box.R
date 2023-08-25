box::use(../ApplicationCode/theme_BMA)
# options('box.path' = getwd())

box::use(
  purrr,                          
  tbl = tibble,                   
  dplyr = dplyr[filter, select],  
  stats[st_filter = filter, ...]  
)

log = function (color) {
  box::use(glue[glue])
  message(glue("My favorite color is {color}"))
}

log('black')

glue('black')
