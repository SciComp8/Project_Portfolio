seurat_obj_2 <- within(seurat_obj@meta.data, {
  log_nFeature_RNA <- log(nFeature_RNA)
  sample_id <- factor(sample_id)
  percent_mt_percentage <- percent.mt * 100
  new_metric <- nCount_RNA / percent.mt  # add a new column 
  rm(orig.ident, nFeature_RNA)  # delete unnecessary columns
})
