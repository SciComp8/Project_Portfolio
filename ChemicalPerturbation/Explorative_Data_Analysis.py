"""
Set up the environment
"""
# pip install pyarrow
# pip install ipython-autotime
# pip install leidenalg
# pip install scanpy

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['figure.figsize'] = [4, 4]
plt.rcParams['figure.dpi'] = 300
sc.set_figure_params(scanpy=True, dpi=100, dpi_save=300, frameon=True, vector_friendly=True, fontsize=14, figsize=None, color_map=None, format='png', facecolor=None, transparent=False, ipython_format='png2x')

"""
Load the preprocessed data
"""
file_path = 'gdrive/MyDrive/Perturbation/'
adata = ad.read_h5ad(file_path + 'single_cell_data/adata_v3_filter.h5ad')
adata

"""
Construct the principal component analysis
"""
sc.tl.pca(adata, svd_solver='arpack')
for gene_i in adata.var_names[:100]:
  sc.pl.pca(adata, color=gene_i)

sc.set_figure_params(scanpy=True, dpi=100, dpi_save=300, frameon=True, vector_friendly=True, fontsize=14, figsize=[6,6], color_map=None, format='png', facecolor=None, transparent=False, ipython_format='png2x')
with plt.rc_context():
  sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
  plt.savefig(file_path + '/analysis_result/PC_number.png', bbox_inches='tight')

adata.write(file_path + 'single_cell_data/adata_pca.h5ad')
adata

"""
Cluster the cells using the neighborhood graph
"""
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.leiden(adata, key_added='leiden_res0_25', resolution=0.25)
sc.tl.leiden(adata, key_added='leiden_res0_5', resolution=0.5)
sc.tl.leiden(adata, key_added='leiden_res1', resolution=1)

sc.set_figure_params(scanpy=True, dpi=100, dpi_save=300, frameon=True, vector_friendly=True, fontsize=14, figsize=[4,4], color_map=None, format='png', facecolor=None, transparent=False, ipython_format='png2x')
for group_i in ['leiden_res0_25', 'leiden_res0_5', 'leiden_res1']:
  print('Resolution: ', group_i)
  sc.tl.paga(adata, groups=group_i)
  with plt.rc_context():
    sc.pl.paga(adata, plot=True, show=False)
    plt.savefig(file_path + '/analysis_result/PAGA_' + group_i + '.png', bbox_inches='tight')

for min_dist_i in [0, 0.25, 0.5, 0.75, 1]:
  print('Minimum distance between embedded points: ', min_dist_i)
  sc.tl.umap(adata, init_pos='paga', min_dist=min_dist_i) 
  sc.pl.umap(adata, color=['leiden_res0_5'], legend_loc='on data', show=False)
  plt.savefig(file_path + '/analysis_result/UMAP_' + str(min_dist_i) + '.png', bbox_inches='tight')

plt.rcParams['figure.figsize'] = [6.4, 4.8]
sc.tl.umap(adata, init_pos='paga', min_dist=0)
sc.pl.umap(adata, color=adata.var_names[:100])

sc.tl.umap(adata, init_pos='paga', min_dist=0.5)
sc.pl.umap(adata, color=adata.var_names[:100], show=True)

adata.write(file_path + 'single_cell_data/adata_cluster.h5ad')
adata

"""
Cluster cells by Cell-type group
"""
sc.tl.paga(adata, groups='cell_type')
plt.rcParams['figure.figsize'] = [4, 4]
with plt.rc_context():
  sc.pl.paga(adata, plot=True, show=False)
  plt.savefig(file_path + '/analysis_result/PAGA_cell_type.png', bbox_inches='tight')

sc.tl.umap(adata, init_pos='paga', min_dist=0)
with plt.rc_context():
  sc.pl.umap(adata, color=['cell_type'], show=False)
  plt.savefig(file_path + '/analysis_result/UMAP_cell_type.png', bbox_inches='tight')

"""
Cluster cells by chemical compound name
"""
sc.tl.paga(adata, groups='sm_name')
with plt.rc_context():
  sc.pl.paga(adata, plot=True, show=False)
  plt.savefig(file_path + '/analysis_result/PAGA_sm_name.png', bbox_inches='tight')

sc.tl.umap(adata, init_pos='paga', min_dist=0)
with plt.rc_context():
  sc.pl.umap(adata, color=['sm_name'], show=False)
  plt.savefig(file_path + '/analysis_result/UMAP_sm_name.png', bbox_inches='tight')

"""
Cluster cells by chemical compound group
"""
