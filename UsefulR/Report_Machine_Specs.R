library(benchmarkme)

get_cpu()$model_name 
get_cpu()$no_of_cores 
get_ram() 
get_platform_info()$OS.type
get_r_version()$platform 
get_r_version()$version.string 

# Print the machine specifications in a file
sink("2024Jan09_machinespecs.txt")
get_cpu()
sink()
