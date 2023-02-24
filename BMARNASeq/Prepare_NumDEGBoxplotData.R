boxplot_data <- function (threshold = 1000, random.seed = 8809678) {
  data.table(threshold=rep(threshold,5), 
             random.seed=random.seed,
             method=c("BMAseq","DESeq2","edgeR","eBayes","voom.limma"),
             variable=rep(c("BMI","AGE","SEX","WBC","BMIxSEX"),each=5),
             num.cDEGs=c(
               intersect(names(BMAseq.eFDR.Main.train$BMI[1:threshold]), names(BMAseq.eFDR.Main.test$BMI[1:threshold])) |> length(),
               intersect(DESeq2.eFDR.GeneName.train$BMI[1:threshold], DESeq2.eFDR.GeneName.test$BMI[1:threshold]) |> length(),
               intersect(edgeR.eFDR.GeneName.train$BMI[1:threshold], edgeR.eFDR.GeneName.test$BMI[1:threshold]) |> length(),
               intersect(names(eBayes.eFDR.train2$BMI[1:threshold]), names(eBayes.eFDR.test2$BMI[1:threshold])) |> length(),
               intersect(names(voom.eFDR.train2$BMI[1:threshold]), names(voom.eFDR.test2$BMI[1:threshold])) |> length(),
               
               intersect(names(BMAseq.eFDR.Main.train$AGE[1:threshold]), names(BMAseq.eFDR.Main.test$AGE[1:threshold])) |> length(),
               intersect(DESeq2.eFDR.GeneName.train$AGE[1:threshold], DESeq2.eFDR.GeneName.test$AGE[1:threshold]) |> length(),
               intersect(edgeR.eFDR.GeneName.train$AGE[1:threshold], edgeR.eFDR.GeneName.test$AGE[1:threshold]) |> length(),
               intersect(names(eBayes.eFDR.train2$AGE[1:threshold]), names(eBayes.eFDR.test2$AGE[1:threshold])) |> length(),
               intersect(names(voom.eFDR.train2$AGE[1:threshold]), names(voom.eFDR.test2$AGE[1:threshold])) |> length(),
               
               intersect(names(BMAseq.eFDR.Main.train$SEX[1:threshold]), names(BMAseq.eFDR.Main.test$SEX[1:threshold])) |> length(),
               intersect(DESeq2.eFDR.GeneName.train$SEX[1:threshold], DESeq2.eFDR.GeneName.test$SEX[1:threshold]) |> length(),
               intersect(edgeR.eFDR.GeneName.train$SEX[1:threshold], edgeR.eFDR.GeneName.test$SEX[1:threshold]) |> length(),
               intersect(names(eBayes.eFDR.train2$SEX[1:threshold]), names(eBayes.eFDR.test2$SEX[1:threshold])) |> length(),
               intersect(names(voom.eFDR.train2$SEX[1:threshold]), names(voom.eFDR.test2$SEX[1:threshold])) |> length(),
               
               intersect(names(BMAseq.eFDR.Main.train$MHABNWBC[1:threshold]), names(BMAseq.eFDR.Main.test$MHABNWBC[1:threshold])) |> length(),
               intersect(DESeq2.eFDR.GeneName.train$MHABNWBC[1:threshold], DESeq2.eFDR.GeneName.test$MHABNWBC[1:threshold]) |> length(),
               intersect(edgeR.eFDR.GeneName.train$MHABNWBC[1:threshold], edgeR.eFDR.GeneName.test$MHABNWBC[1:threshold]) |> length(),
               intersect(names(eBayes.eFDR.train2$MHABNWBC[1:threshold]), names(eBayes.eFDR.test2$MHABNWBC[1:threshold])) |> length(),
               intersect(names(voom.eFDR.train2$MHABNWBC[1:threshold]), names(voom.eFDR.test2$MHABNWBC[1:threshold])) |> length(),
               
               intersect(names(BMAseq.eFDR.Interaction.train$BMI[1:threshold]), names(BMAseq.eFDR.Interaction.test$BMI[1:threshold])) |> length(),
               intersect(DESeq2.eFDR.GeneName.train$BMIxSEX[1:threshold], DESeq2.eFDR.GeneName.test$BMIxSEX[1:threshold]) |> length(),
               intersect(edgeR.eFDR.GeneName.train$BMIxSEX[1:threshold], edgeR.eFDR.GeneName.test$BMIxSEX[1:threshold]) |> length(),
               intersect(names(eBayes.eFDR.train2$BMIxSEX[1:threshold]), names(eBayes.eFDR.test2$BMIxSEX[1:threshold])) |> length(),
               intersect(names(voom.eFDR.train2$BMIxSEX[1:threshold]), names(voom.eFDR.test2$BMIxSEX[1:threshold])) |> length()))
}
