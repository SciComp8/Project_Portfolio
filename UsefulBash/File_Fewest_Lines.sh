# Find the file with the fewest lines
wc -l relative_path/*.csv | grep -v total | sort -n | head -n 1
