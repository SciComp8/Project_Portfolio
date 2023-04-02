# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
# Ref: https://cran.r-project.org/web/packages/xfun/vignettes/xfun.html

# Install the self-made library
# devtools::install("/Library/Frameworks/R.framework/Versions/4.2/Resources/library/knitautomator")


knit_filename <- function(input, ...) {
  out_dir <- "../report/";
  temp <- sub(".*\\d_", "", xfun::sans_ext(input));
  temp <- sub("_.*", "", temp); 
  out_project <- paste0("_", temp, "_"); 
  out_author <- "ALiu";
  if (grepl("ex", xfun::sans_ext(input)) == T) {
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
