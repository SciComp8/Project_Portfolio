## Load all data
data_patient <- read.delim("../data/raw/data_clinical_patient.txt")
data_patient <- data_patient[-c(1:4), ]
names(data_patient) <- 
  sapply(strsplit(names(data_patient), split = "\\."), function(x)
    paste(x[x != ""], collapse = "."))

data_sample <- read.delim("../data/raw/data_clinical_sample.txt")
data_sample <- data_sample[-c(1:4), ]
names(data_sample) <- 
  sapply(strsplit(names(data_sample), split = "\\."), function(x)
    paste(x[x != ""], collapse = "."))

data_genepanel <- read.delim("../data/raw/data_gene_panel_matrix.txt") # not valuable

data_mutation <- read.delim("../data/raw/data_mutations.txt") # tumor sample barbode = sample identifier 

data_sv <- read.delim("../data/raw/data_sv.txt") # sample identifier
# unique(data_sv$Sample_Id) # 285 unique samples 

## Merge all necessary data
data_full <- merge(x = data_patient, 
                   y = data_sample, 
                   by.x = "X.Patient.Identifier", 
                   by.y = "X.Patient.Identifier")

data_full <- merge(x = data_full,
                   y = data_mutation,
                   by.x = "Sample.Identifier",
                   by.y = "Tumor_Sample_Barcode")

data_full <- merge(x = data_full,
                   y = data_sv,
                   by.x = "Sample.Identifier",
                   by.y = "Sample_Id")

## Write the derived data
dput(data_full, "../data/derived/2022Nov8_1661tmb.Rdata")
write.csv(data_full, "../data/derived/2022Nov8_1661tmb.csv", row.names = F)
