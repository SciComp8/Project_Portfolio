import numpy as np 
import pandas as pd # pip install pyarrow
import scanpy as sc # pip install numpy==1.25
import anndata as ad
from scipy.sparse import csr_matrix
import scipy.sparse
from scipy.stats import median_abs_deviation

"""
Load the raw count data
"""
file_path = '../single_cell_data/'
adata_train_raw = pd.read_parquet(file_path + 'adata_train.parquet', engine='pyarrow')
print(adata_train_raw.shape)

"""
Transform and compress the raw count data to the cell x gene sparse matrix
"""
adata_train_raw_wide = adata_train_raw.pivot(index='obs_id', columns='gene', values='count')
adata_train_raw_wide.head()
adata_train_raw_wide.to_parquet(file_path + 'count_matrix_train.parquet')
count_matrix = pd.read_parquet(file_path + 'count_matrix_train.parquet', engine='pyarrow') 
count_matrix.head()
count_matrix_new = count_matrix.fillna(0)
count_matrix_new.head()
col_label = count_matrix_new.columns 
count_matrix_new.to_parquet(file_path + 'count_matrix_train_new.parquet', engine='pyarrow')
count_matrix_sparse = count_matrix_new.astype(pd.SparseDtype("float64",0))
count_matrix_csr = count_matrix_sparse.sparse.to_coo().tocsr()
scipy.sparse.save_npz(file_path + 'count_matrix_csr_train.npz', count_matrix_csr)
sparse_count_matrix = scipy.sparse.load_npz(file_path + 'count_matrix_csr_train.npz') 
print(sparse_count_matrix.shape)

"""
Load the meta data
"""
df_meta_train = pd.read_csv(file_path + 'adata_obs_meta_new.csv')
for i in [3, 5, 7, 9]:
    df_meta_train['sm_cluster_' + str(i)] = pd.Categorical(df_meta_train['sm_cluster_' + str(i)])
print(df_meta_train.shape) 
df_meta_train.head()
df_meta_train = df_meta_train.set_index('obs_id')

"""
Build an AnnData object
"""
adata = sc.AnnData(X=sparse_count_matrix, obs=df_meta_train)
adata

adata.var_names = col_label
adata.write(file_path + 'adata_v1.h5ad')

"""
Perform the cell quality control
"""
adata = ad.read_h5ad(file_path + 'adata_v1.h5ad')
adata[0:6,0:6].to_df()
# adata.X

sc.pl.highest_expr_genes(adata, n_top=10)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata

adata.var['mt'] = adata.var_names.str.startswith('MT-')
np.count_nonzero(adata.var['mt']) # 37
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[20], log1p=True, inplace=True)
adata

sc.pl.violin(adata, keys='n_genes_by_counts', jitter=0.5)
sc.pl.violin(adata, keys='total_counts', jitter=0.5)
sc.pl.violin(adata, keys='pct_counts_mt', jitter=0.5)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)

adata.obs.outlier.value_counts()

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)

adata.obs.mt_outlier.value_counts()

print(f"Total number of cells: {adata.n_obs}")
adata_2 = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Number of cells after filtering of low quality cells: {adata_2.n_obs}") 

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

adata_2.write(file_path + 'adata_v2.h5ad')

"""
Normalize the count matrix
"""
adata = ad.read_h5ad(file_path + 'adata_v2.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

"""
Select the highly variable genes (feature selection)
"""
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.var.highly_variable
adata.var.highly_variable.value_counts()
adata.raw = adata
# adata_back = adata.raw.to_adata()
adata = adata[:, adata.var.highly_variable]
adata

"""
Regress out effects of total counts and percentage of mitochondrial genes, and scale the data to the unit variance
"""
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
adata[0:6,0:6].to_df()

"""
Save the final preprocessed file
"""
adata.write(file_path + 'adata_v3_filter.h5ad')
