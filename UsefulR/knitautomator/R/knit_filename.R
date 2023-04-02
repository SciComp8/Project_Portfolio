# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

knit_filename <- function(input, ...) {
  out_dir <- "../report/";
  temp <- sub(".*\\d_", "", xfun::sans_ext(input)); # Newly add
  temp <- sub("_.*", "", temp); # Newly add
  out_project <- paste0("_", temp, "_"); # Newly add
  out_author <- "ALiu_XKZ";
  if (grepl("ex", xfun::sans_ext(input)) == T) { # Newly add
    out_v <- "_V1";
  } else {
    out_v <- "_V2";
  }
  file_type <- ".docx";
  rmarkdown::render(input,
                    encoding = encoding,
                    output_file = paste0(
                      out_dir,
                      format(Sys.Date(), "%Y%b%d"),
                      out_project,
                      out_author,
                      out_v,
                      file_type),
                    envir = globalenv())
}
