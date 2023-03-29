library(optparse)

# Define a list of options using `make_option()`
option_list = list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = "File path to be processed"))

# Create an `OptionParser` object and pass in the option list 
opt_parser = OptionParser(option_list = option_list)

# Use `parse_args()` to parse the command-line arguments
opt = parse_args(opt_parser)

# Extract the option values
file_path = opt$file

# Check if the file option was specified, and if so, we do something with the file.
if (is.null(file_path)) {
  stop("No file path specified.")
} else {
  print(file_path)
}

# Test optparse in the command line
# anniliu@Angela-Excelsior-MacBook Desktop % Rscript optparse.R -f "C:\\D"
# [1] "C:\\D"
