# pip install rdkit
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.cluster import AgglomerativeClustering
import numpy as np

"""
Load the meta dataset into a pandas DataFrame
"""
file_path = '../single_cell_data/'
df = pd.read_csv(file_path + 'adata_obs_meta.csv')
df_SMILES_unique = df["SMILES"].unique()
print(df_SMILES_unique[:6])

"""
Define a function to calculate Tanimoto similarity between two SMILES strings
"""
def calculate_similarity(smiles1, smiles2, radius=2, nBits=2048):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return None
    fingerprint1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=nBits)
    fingerprint2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=nBits)
    similarity = DataStructs.TanimotoSimilarity(fingerprint1, fingerprint2)
    return similarity

"""
Calculate pairwise similarity matrix
"""
num_compounds = len(df_SMILES_unique)
similarity_matrix = np.zeros((num_compounds, num_compounds))

for i in range(num_compounds):
    for j in range(i, num_compounds):
        smiles1 = df_SMILES_unique[i]
        smiles2 = df_SMILES_unique[j]
        similarity = calculate_similarity(smiles1, smiles2)
        similarity_matrix[i, j] = similarity
        similarity_matrix[j, i] = similarity # The returned similarity matrix is symmetric and the elements are in [0, 1]

distance_matrix = 1 - similarity_matrix

"""
Cluster compounds into k clusters using Agglomerative Clustering
"""
cluster_label_list = [AgglomerativeClustering(n_clusters=i, linkage="ward").fit_predict(distance_matrix) for i in [3, 5, 7, 9]]

"""
Add cluster labels to the DataFrame
"""
cluster_label_list
result = [dict(zip(df_SMILES_unique, cluster_label_list[i])) for i in range(4)]
df['sm_cluster_3'] = df['SMILES'].map(result[0])
df['sm_cluster_5'] = df['SMILES'].map(result[1])
df['sm_cluster_7'] = df['SMILES'].map(result[2])
df['sm_cluster_9'] = df['SMILES'].map(result[3])

"""
Veiw the distribution of clustering labels
"""
[df['sm_cluster_' + str(i)].value_counts() for i in [3, 5, 7, 9]]

"""
Save the DataFrame with the new compound labels
"""
df.to_csv(file_path + "adata_obs_meta_new.csv", index=False)
