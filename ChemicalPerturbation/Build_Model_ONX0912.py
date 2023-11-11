"""
Predict the perturbation effect of onx0912 
"""

pip install anndata==0.9.2
pip install git+https://github.com/theislab/scgen.git
pip install pertpy
pip install scanpy
pip install numba
pip install llvmlite

import scanpy as sc
import pertpy as pt
import scgen
import anndata as ad
import numpy as np

# Load the data
from google.colab import drive
drive.mount('/content/gdrive')

file_path = 'gdrive/MyDrive/Perturbation/'
adata = ad.read_h5ad(file_path + 'single_cell_data/adata_v3_filter.h5ad')
adata

# Count per categorical level
def count_entries(df, *args):
    """Return a dictionary with counts of occurrences as value for each key."""
    cols_count = {}
    for col_name in args:
        col = df[col_name]
        for entry in col:
            if entry in cols_count.keys():
                cols_count[entry] += 1
            else:
                cols_count[entry] = 1
    return cols_count

df = adata.obs
sm_name_count = count_entries(df, 'sm_name')
print(sm_name_count)
sm_cluster_5_count = count_entries(df, 'sm_cluster_5')
print(sm_cluster_5_count)

# ONX 0912
df_col = adata.obs['sm_name']
is_ONX0912 = df_col.str.contains('ONX 0912')
df_col[is_ONX0912]
# Format the condition types
adata.obs['condition'] = ['perturbated' if sm_i == 'Oprozomib (ONX 0912)' else 'control' for sm_i in adata.obs['sm_name']]

# Count the condition types
df = adata.obs
condition_count = count_entries(df, 'condition')
print(condition_count)

# Count the cell types
adata.obs.cell_type.value_counts()

# Remove all perturbated CD4+T and CD8+T cells from the training data
adata_train = adata[
    ~(
        (adata.obs['cell_type'] == 'T cells CD4+')
        & (adata.obs['condition'] == 'perturbated')
    )
].copy()
adata_train.obs.cell_type.value_counts()

# Keep all perturbated CD4+T and CD8+T cells as the training set
cd4_perturb = adata[
    (
        (adata.obs['cell_type'] == 'T cells CD4+')
        & (adata.obs['condition'] == 'perturbated')
    )
].copy()
cd4_perturb.obs.cell_type.value_counts()

scgen.SCGEN.setup_anndata(adata_train, batch_key='condition', labels_key='cell_type')

# Save the data
adata_train.write(file_path + 'single_cell_data/adata_train.h5ad')
cd4_perturb.write(file_path + 'single_cell_data/cd4_perturb.h5ad')

# Train the model
onx0912_model = scgen.SCGEN(adata_train, n_hidden=800, n_latent=100, n_layers=2)
onx0912_model.train(
    max_epochs=100, batch_size=32, early_stopping=True, early_stopping_patience=25
)
# INFO:pytorch_lightning.utilities.rank_zero:GPU available: True (cuda), used: True
# INFO:pytorch_lightning.utilities.rank_zero:TPU available: False, using: 0 TPU cores
# INFO:pytorch_lightning.utilities.rank_zero:IPU available: False, using: 0 IPUs
# INFO:pytorch_lightning.utilities.rank_zero:HPU available: False, using: 0 HPUs
# INFO: LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]
# INFO:lightning.pytorch.accelerators.cuda:LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]
# Epoch 26/100:  26%|██▌       | 26/100 [06:38<18:54, 15.33s/it, v_num=1, train_loss_step=615, train_loss_epoch=468]
# Monitored metric elbo_validation did not improve in the last 25 records. Best score: 1902.328. Signaling Trainer to stop.

# Get the 100-dimensional vector for each cell
adata_train.obsm['scgen'] = onx0912_model.get_latent_representation()

# Visualize the new embedding from the latent representation in a UMAP plot
sc.pp.neighbors(adata_train, use_rep='scgen')
sc.tl.umap(adata_train)
sc.pl.umap(adata_train, color=['condition', 'cell_type'], wspace=0.4, frameon=False)

# Predict the CD4+T responses to ONX 0912
pred, delta = onx0912_model.predict(
    ctrl_key='control', stim_key='perturbated', celltype_to_predict='T cells CD4+'
)

pred.obs['condition'] = 'predicted perturbated'
