import scanpy as sc
import pertpy as pt
import scgen
import anndata as ad
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from adjustText import adjust_text
from google.colab import drive

class DrugEffectPredictor:
    def __init__(self, file_path, drug_name):
        self.file_path = file_path
        self.drug_name = drug_name
        self.adata = None
        self.adata_train = None
        self.cd4_perturb = None
        self.model = None
        self.adata_evaluate = None
        self.deg = None

    def load_data(self):
        drive.mount('/content/gdrive')
        self.adata = ad.read_h5ad(self.file_path + 'single_cell_data/adata_v3_filter.h5ad')

    def preprocess_data(self):
        self.adata.obs['condition'] = ['perturbated' if self.drug_name in sm_name else 'control' for sm_name in self.adata.obs['sm_name']]
        self.adata_train = self.adata[~((self.adata.obs['cell_type'] == 'T cells CD4+') & (self.adata.obs['condition'] == 'perturbated'))].copy()
        self.cd4_perturb = self.adata[((self.adata.obs['cell_type'] == 'T cells CD4+') & (self.adata.obs['condition'] == 'perturbated'))].copy()

    def train_model(self):
        scgen.SCGEN.setup_anndata(self.adata_train, batch_key='condition', labels_key='cell_type')
        self.model = scgen.SCGEN(self.adata_train, n_hidden=800, n_latent=100, n_layers=2)
        self.model.train(max_epochs=100, batch_size=32, early_stopping=True, early_stopping_patience=25)
        self.adata_train.obsm['scgen'] = self.model.get_latent_representation()
        self.adata_train.write(self.file_path + 'single_cell_data/adata_train.h5ad')

    def predict_and_evaluate(self):
        pred, delta = self.model.predict(ctrl_key='control', stim_key='perturbated', celltype_to_predict='T cells CD4+')
        pred.obs['condition'] = 'predicted perturbated'
        adata_control = self.adata[((self.adata.obs['cell_type'] == 'T cells CD4+') & (self.adata.obs['condition'] == 'control'))]
        self.adata_evaluate = adata_control.concatenate(self.cd4_perturb, pred)
        sc.tl.pca(self.adata_evaluate)
        sc.pl.pca(self.adata_evaluate, color='condition', frameon=False)
        self.adata_cd4 = self.adata[self.adata.obs['cell_type'] == 'T cells CD4+']
        sc.tl.rank_genes_groups(self.adata_cd4, groupby='condition', method='wilcoxon')
        self.deg = self.adata_cd4.uns['rank_genes_groups']['names']['perturbated']

    def visualize_results(self):
        deg_list = np.intersect1d(deg, adata_evaluate.var_names)
    
        # Calculate the R² correlation between mean gene expression of predicted and existing CD4+T cells
        pred_deg = pred[:, deg_list]
        true_deg = cd4_perturb[:, deg_list]
        x_deg = np.asarray(np.mean(pred_deg.X, axis=0)).ravel()
        y_deg = np.asarray(np.mean(true_deg.X, axis=0)).ravel()
        slope, intercept, r_value_deg, p_value_deg, std_err_deg = stats.linregress(x_deg, y_deg)
        r_squared = r_value_deg**2
    
        # Visualize the R² correlation
        plt.figure(figsize=(10, 8))
        plt.scatter(x_deg, y_deg, alpha=0.7)
        plt.plot(x_deg, intercept + slope * x_deg, 'r', label=f'y = {slope:.2f}x + {intercept:.2f}')
        texts = [plt.text(x_deg[i], y_deg[i], deg_list[i], fontsize=9) for i in range(min(len(deg_list), 10))]
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), expand_text=(2, 3))
        plt.text(min(x_deg), max(y_deg), f'$R^2 = {r_squared:.3f}$', fontsize=12)
        plt.xlabel('Predicted mean gene expression in CD4+T control cells')
        plt.ylabel('True mean gene expression in CD4+T perturbated cells')
        plt.title('Scatter plot of the mean gene expression between prediction and ground truth')
        # plt.legend()
        plt.show()
    
        # Visualize the distribution of the top differentially expressed genes
        for gene in deg_list[:10]:  # Adjust the number of genes as needed
            sc.pl.violin(adata_evaluate, gene, groupby='condition')

# Use the class
file_path = 'gdrive/MyDrive/Perturbation/'
drug_name = 'K-02288' # Replace with the drug we want to analyze

predictor = DrugEffectPredictor(file_path, drug_name)
predictor.load_data()
predictor.preprocess_data()
predictor.train_model()
predictor.predict_and_evaluate()
predictor.visualize_results()

