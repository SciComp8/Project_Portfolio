dat.work <- dat0[dat0$SxRT == "LumpectomyRTYes" | dat0$SxRT == "LumpectomyRTNo", ] # Not use which() here; but mis-classify 155 patients in 0 or 1 recur group into NA group

table(SxRT = dat0$SxRT, Recur = dat0$Recur, useNA = "ifany")
#                  Recur
# SxRT                0   1
#   LumpectomyRTYes 454  36
#   OtherSxRTYes     30   0
#   LumpectomyRTNo  120  21
#   OtherSxRTNo     180   7
#   <NA>            151   4

table(SxRT = dat.work$SxRT, Recur = dat.work$Recur, useNA = "ifany")
#                  Recur
# SxRT                0   1 <NA>
#   LumpectomyRTYes 454  36    0
#   OtherSxRTYes      0   0    0
#   LumpectomyRTNo  120  21    0
#   OtherSxRTNo       0   0    0
#   <NA>              0   0  155
  
non.miss.idx <- which(dat0$SxRT == "LumpectomyRTYes" | dat0$SxRT == "LumpectomyRTNo")
miss.idx <- which(is.na(dat0$SxRT))
dat.work <- dat0[c(non.miss.idx, miss.idx), ]
table(SxRT = dat.work$SxRT, Recur = dat.work$Recur, useNA = "ifany")
#                  Recur
# SxRT                0   1
#   LumpectomyRTYes 454  36
#   OtherSxRTYes      0   0
#   LumpectomyRTNo  120  21
#   OtherSxRTNo       0   0
#   <NA>            151   4
  
dat.work$SxRT <- droplevels(dat.work$SxRT)
table(dat.work$SxRT, useNA = "ifany")
# LumpectomyRTYes  LumpectomyRTNo            <NA> 
#             490             141             155 

DDX1.pos.train <- which(grepl(pattern = "ENSG00000079785", x = names(temp$BMAseq.eFDR.Main.train[["BMI"]])))
