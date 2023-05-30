names(data_full)
 [1] "Patient ID"                              
 [2] "Patient - Gender"                        
 [3] "Age at diagnosis"                        
 [4] "Race/Ethnicity"                          
 [5] "Date of diagnosis"                       
 [6] "Procedure"                               
 [7] "Axila Procedure"                         
 [8] "Chemotherapy"                            
 [9] "Endocrine therapy"                       
[10] "Radiation therapy"                       
[11] "BMI at diagnosis"                        
[12] "Height (m)"                              
[13] "Weight at diagnosis (kg)"                
[14] "Pathologic T stage"                      
[15] "Tumor size (cm)"                         
[16] "N stage"                                 
[17] "ER Status"                               
[18] "PR Status"                               
[19] "HER2 Status"                             
[20] "Oncotype Score"                          
[21] "Histological type"                       
[22] "Histological grade"                      
[23] "Lymphovascular invasion"                 
[24] "Date of last follow-up"                  
[25] "Local recurrence at follow-up"           
[26] "Date of local recurrence (if applicable)"
[27] "Distant metastasis at follow-up"         
[28] "Date of distant metastasis"              
[29] "Duration of follow-up in days"     

##------Format variable names into the unified style (e.g., Patient ID -> Patient.id; Race/Ethnicity -> Race.ethnicity)------
temp_name <- names(data_full)
temp_name <- sapply(strsplit(temp_name, split = "\\-|\\s|\\(|\\)|\\/|\\%|^the$"), function(x) paste(x[!x %in% c("", "-", "(", ")", "/", "%", "the")], collapse = ".")) # \\s means the literal space 
temp_name <- gsub("$", ".", temp_name) # Chemotherapy without any comma [differ from Patient.ID]
temp_name <- sub("\\..*", ".", temp_name) %0% tolower(substring(temp_name, regexpr("*.\\.", temp_name) + 2))
temp_name <- gsub(".$", "", temp_name)
names(data_full) <- temp_name

