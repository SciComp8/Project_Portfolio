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

# Load the data
from google.colab import drive
drive.mount('/content/gdrive')

file_path = 'gdrive/MyDrive/Perturbation/'
adata = ad.read_h5ad(file_path + 'single_cell_data/adata_v3_filter.h5ad')
adata

def count_entries(df, *args):
    """Return a dictionary with counts of occurrences as value for each key."""

    # Initialize an empty dictionary: cols_count
    cols_count = {}

    # Iterate over column names in args
    for col_name in args:

        # Extract column from DataFrame: col
        col = df[col_name]

        # Iterate over the column in DataFrame
        for entry in col:

            # If entry is in cols_count, add 1
            if entry in cols_count.keys():
                cols_count[entry] += 1

            # Else add the entry to cols_count, set the value to 1
            else:
                cols_count[entry] = 1

    # Return the cols_count dictionary
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
        (adata.obs[adata.obs["cell_type"].isin(["T cells CD4+", "T cells CD8+"])])
        & (adata.obs["condition"] == "perturbated")
    )
].copy()
adata_train.obs.cell_type.value_counts()

# Keep all perturbated CD4+T and CD8+T cells as the training set
cd48_perturb = adata[
    (
        (adata.obs[adata.obs["cell_type"].isin(["T cells CD4+", "T cells CD8+"])])
        & (adata.obs["condition"] == "perturbated")
    )
].copy()
cd48_perturb.obs.cell_type.value_counts()

scgen.SCGEN.setup_anndata(adata_train, batch_key="condition", labels_key="cell_type")

# Save the data
adata_train.write(file_path + 'single_cell_data/adata_train.h5ad')
cd48_perturb.write(file_path + 'single_cell_data/cd48_perturb.h5ad')
