# Find the file with the fewest lines
wc -l seasonal/*.csv | grep -v total | sort -n | head -n 1
