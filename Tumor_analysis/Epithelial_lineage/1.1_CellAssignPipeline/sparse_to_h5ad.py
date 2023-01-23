import anndata
from scipy import sparse, io
import scanpy as sc
import pandas as pd
Import numpy as np

meta_data = pd.read_csv("IntegratedTumorMetaData.csv")

matrix = io.mmread('IntegratedTumorsRNARaw.mtx')
matrix = matrix.todense()
integrated_adata = anndata.AnnData(X=matrix)
integrated_adata = integrated_adata.T
integrated_adata.obs = meta_data
   
sc.pp.normalize_total(integrated_adata, target_sum=1e4)

lib_size = integrated_adata.X.sum(1)

integrated_adata.obs["size_factor"] = lib_size / np.mean(lib_size)

integrated_adata = integrated_adata.copy()


gene_names = pd.read_csv('gene_names.csv', header=None)
integrated_adata.var = gene_names

sigs_excel = "Epithelial_lung_signature.xlsx"
covid_sigs_df = pd.read_excel(sigs_excel).T
covid_sigs = {}
for i,row in covid_sigs_df.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(sig.strip("\xa0").strip())
    covid_sigs[i] = list(set(sig_set))
covid_sigs

laughney_sigs_df = pd.read_excel(sigs_excel, sheet_name='Laughney_signature').T
laughney_sigs = {}
for i,row in laughney_sigs_df.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(str(sig).strip("\xa0").strip())
    laughney_sigs[i] = list(set(sig_set))
laughney_sigs


epi_sigs_df = pd.read_excel(sigs_excel, sheet_name='Epithelial_lung').T
epi_sigs = {}
for i,row in epi_sigs_df.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(str(sig).strip("\xa0").strip())
    epi_sigs[i] = list(set(sig_set))

all_sigs = {
    "covid": covid_sigs,
    "laughney": laughney_sigs,
    "epi": epi_sigs
}

def create_marker_gene_matrix(sigs_dict):
    all_genes = []
    for x in sigs_dict:
        all_genes.extend(sigs_dict[x])
    all_genes = list(set(all_genes))
    df = []
    for x in sigs_dict:
        df_col = []
        for y in all_genes:
            df_col.append(int( y in sigs_dict[x]))
        df.append(df_col)
    return pd.DataFrame(df, index =list(sigs_dict.keys()), columns = all_genes).T

gene_mat1 = create_marker_gene_matrix(covid_sigs)
gene_mat2 = create_marker_gene_matrix(laughney_sigs)
gene_mat3 = create_marker_gene_matrix(covid_sigs)

all_genes = set(gene_mat1.index)& set(gene_mat2.index)&set(gene_mat3.index)

integrated_adata.var["gene-name"]= integrated_adata.var[0]
integrated_adata.var = integrated_adata.var.set_index(0)

integrated_adata.var = integrated_adata.var.set_index("gene-name")
integrated_adata_filt = integrated_adata[:,integrated_adata.var.index.isin(all_genes)]
integrated_adata_filt.var.rename(columns={0:"genes"}, inplace= True)


integrated_adata_filt.write(“filtered_cellassign_adata.h5ad”)

integrated_adata.write("full_integrated_adata.h5ad")
