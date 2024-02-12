## Example usage 1
library(optparse)

# Define a list of options using `make_option()`
option_list = list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = "File path to be processed"))

# Create an `OptionParser` object and pass in the option list 
opt_parser = OptionParser(option_list = option_list)

# Use `parse_args()` to parse the command-line arguments
opt = parse_args(opt_parser)

# ! Ensure optparse is correctly parsing arguments
if(exists("opt")) {
  print("opt exists")
} else {
  print("opt does not exist")
}

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


## Example usage 2
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
  # do something with the file
  dat0 <- fst::read_fst(path = file_path)
  all.char <- sapply(names(dat0)[sapply(dat0, is.character)], function(x) with(dat0, table(get(x), useNA = "ifany"))) 
  all.fac <- sapply(names(dat0)[sapply(dat0, is.factor)], function(x) with(dat0, table(get(x), useNA = "ifany"))) 
  all.cat <- list(all.char, all.fac)
  all.cat
}

# Test optparse in the terminal 2
# Rscript optparse.R -f /Users/Documents/data/derived/TNBC.RData

                    
