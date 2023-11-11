from scipy import stats
import matplotlib.pyplot as plt
from adjustText import adjust_text

adata_control = adata[
    ((adata.obs['cell_type'] == 'T cells CD4+') & (adata.obs['condition'] == 'control'))
]

# Combine control CD4+T cell, predicted CD4+T cell, truly perturbated CD4+T cells
adata_evaluate = adata_control.concatenate(cd4_perturb, pred)

# Build the PCA to embed and compare the latent spaces of gene expression in CD4+T control cells (Tc), CD4+T truly perturbed cells (Tt), CD4+T predicted cells (Tp)
sc.tl.pca(adata_evaluate)
sc.pl.pca(adata_evaluate, color='condition', frameon=False)

adata_cd4 = adata[adata.obs['cell_type'] == 'T cells CD4+']

# Detect the differentially expressed genes
sc.tl.rank_genes_groups(adata_cd4, groupby='condition', method='wilcoxon')
deg = adata_cd4.uns['rank_genes_groups']['names']['perturbated']
deg

deg_list = np.intersect1d(deg, adata_evaluate.var_names)

# Calculate the R² correlation between mean gene expression of predicted and existing CD4+T cells
pred_deg = pred[:, deg_list]
true_deg = cd4_perturb[:, deg_list]
x_deg = np.asarray(np.mean(pred_deg.X, axis=0)).ravel()
y_deg = np.asarray(np.mean(true_deg.X, axis=0)).ravel()
slope, intercept, r_value_deg, p_value_deg, std_err_deg = stats.linregress(
  x_deg, y_deg)
r_squared = r_value_deg**2

# Visualize the R² correlation between mean gene expression of predicted and existing CD4+T cells
plt.figure(figsize=(10, 8))
plt.scatter(x_deg, y_deg, alpha=0.7)
for idx in range(10):
    plt.annotate(deg_list[idx], (x_deg[idx], y_deg[idx]), fontsize=12)
plt.plot(x_deg, intercept + slope * x_deg, 'r', label=f'y = {slope:.2f}x + {intercept:.2f}')
texts = []
for idx in range(10):
    texts.append(plt.text(x_deg[idx], y_deg[idx], deg_list[idx], fontsize=9))

adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), expand_text=(2, 3))
plt.text(min(x_deg), max(y_deg), f'$R^2 = {r_squared:.3f}$', fontsize=12)
plt.xlabel('Predicted mean gene expression in the CD4+T control cells')
plt.ylabel('True mean gene expression in the CD4+T perturbated cells')
plt.title('Scatter plot of the mean gene expression between prediction and ground truth')
plt.show()

# Visualize the distribution of the top differentially expressed genes between control CD4+T cells and condition CD4+T cells
for gene in deg_list[:10]:
    sc.pl.violin(adata_evaluate, gene, groupby='condition')

