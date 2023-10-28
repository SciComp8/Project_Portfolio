paste -d "|" file_name_1 file_name_2 file_name_3

paste -d "|," file_name_1 file_name_2 file_name_3
# The contents between the 1st and 2nd files are separated by |
# The contents between the 2nd and 3rd files are separated by ,

paste -s -d ":" file_name_1 file_name_2 file_name_3
# -s: read all the lines from a single file and merge all these lines into a single line with each line separated by tab; these single lines are separated by newline

cut -d " " -f 1 file_name_2 | paste file_name_1 -
