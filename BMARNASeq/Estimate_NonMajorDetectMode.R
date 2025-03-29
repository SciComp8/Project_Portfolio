BMI_5000_ClassFreq <- readRDS(file='../../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/BMI_5000_ClassFreq')
cDEG_Num_NoMajorConsist <- filter(.data=BMI_5000_ClassFreq, !Class%in%c('100000000', '010000000', '000100000', '010100000', '001010000', '000110000', '100000101', '011110000', '100001111', '101011111', '111111111')) |>
  group_by(Seed)
write.csv(cDEG_Num_NoMajorConsist, '../../ApplicationData/derived/RandomSeed/HeatmapBoxplotData/BMI_cDEG_Num_NoMajorConsist.csv')

cDEG_Num_NoMajorConsist_BySeed <- filter(.data=BMI_5000_ClassFreq, !Class%in%c('100000000', '010000000', '000100000', '010100000', '001010000', '000110000', '100000101', '011110000', '100001111', '101011111', '111111111')) |>
  group_by(Seed) |>
  summarise(cDEG_Num = sum(Frequency))

summary(cDEG_Num_NoMajorConsist_BySeed$cDEG_Num)[3]
summary(cDEG_Num_NoMajorConsist_BySeed$cDEG_Num)[1]
summary(cDEG_Num_NoMajorConsist_BySeed$cDEG_Num)[6]
