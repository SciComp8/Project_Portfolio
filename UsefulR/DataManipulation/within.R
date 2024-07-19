seurat_obj_2 <- within(seurat_obj@meta.data, {
  log_nFeature_RNA <- log(nFeature_RNA) # log-transform a metadata column
  sample_id <- factor(sample_id) # covert a numeric metadata column to a factor
  percent_mt_percentage <- percent.mt * 100 # add a new column 
  new_metric <- nCount_RNA / percent.mt  # add a new column 
  rm(orig.ident, nFeature_RNA)  # delete unnecessary columns
})

# within() examines the environment after the evaluation of expr 
# and makes the corresponding modifications to a copy of data (this may fail in the data frame case if objects are created which cannot be stored in a data frame), and returns it. 
# within() can be used as an alternative to transform.
