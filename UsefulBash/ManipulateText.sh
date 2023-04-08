data=$1
cat  $data | grep non-missing | awk '{if ($2=="Asian") print $6, $7}' | sed 's/"//g' | sed 's/://g' |sort -u > newdata.csv

# grep non-missing: filters out any lines that do not contain the string non-missing.
# awk '{if ($2=="Asian") print $6, $7}': selects the 6th and 7th fields of each line if the 2nd field is "Asian".
# sed 's/"//g': removes all double quotes from the selected fields.
# sed 's/://g' removes all colons from the selected fields.
# sort -u: sorts the resulting two fields alphabetically and removes duplicates.
