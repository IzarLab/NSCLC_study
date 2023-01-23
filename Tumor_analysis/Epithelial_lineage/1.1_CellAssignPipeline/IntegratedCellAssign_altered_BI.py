#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import numpy as np
import scvi
import logging
import subprocess


integrated_save_file = 'data/IntegratedDataQCCheckpoints/integrated_normalized_control.h5ad'



# In[2]:


#Read All Signatures, Store as Dictionary
covid_sigs_df = pd.read_excel("Epithelial_lung_signature_BI_ST.xlsx").T
covid_sigs = {}
for i,row in covid_sigs_df.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(sig.strip("\xa0").strip())
    covid_sigs[i] = list(set(sig_set))
covid_sigs


#Read All Signatures, Store as Dictionary
covid_sigs_df2 = pd.read_excel("Epithelial_lung_signature_BI.xlsx").T
covid_sigs2 = {}
for i,row in covid_sigs_df2.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(sig.strip("\xa0").strip())
    covid_sigs2[i] = list(set(sig_set))
covid_sigs2


laughney_sigs_df = pd.read_excel("Epithelial_lung_signature_BI.xlsx", sheet_name='Developmental signatures').T
laughney_sigs = {}
for i,row in laughney_sigs_df.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(str(sig).strip("\xa0").strip())
    laughney_sigs[i] = list(set(sig_set))
laughney_sigs


epi_sigs_df = pd.read_excel("Epithelial_lung_signature_BI.xlsx", sheet_name='CellAssign_signatures long ').T
epi_sigs = {}
for i,row in epi_sigs_df.fillna(0).iterrows():
    sig_set = []
    for sig in row:
        if sig != 0:
            sig_set.append(str(sig).strip("\xa0").strip())
    epi_sigs[i] = list(set(sig_set))

all_sigs = {
    "dev": covid_sigs,
    #"dev": laughney_sigs,
    #"long": epi_sigs
}

def create_marker_gene_matrix(sigs_dict , sample_genes = []):
    if (len(sample_genes) == 0):
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
    else:
        all_genes = []
        for x in sigs_dict:
            all_genes.extend(sigs_dict[x])
        all_genes = list(set(all_genes) & set(sample_genes))
        df = []
        for x in sigs_dict:
            df_col = []
            for y in all_genes:
                df_col.append(int( y in sigs_dict[x]))
            df.append(df_col)
        return pd.DataFrame(df, index =list(sigs_dict.keys()), columns = all_genes).T

create_marker_gene_matrix(covid_sigs)


def get_cellassign_sample(sample, matrix):
    gene_matrix =  create_marker_gene_matrix(matrix)
    sample= sample[:,sample.var.index.isin(set(gene_matrix.index))]
    sample = sample.copy()
    return sample


# In[3]:


adata = sc.read_h5ad("filtered_cellassign_adata_ST.h5ad")



adata_filts = {}
for y in list(all_sigs.keys()):
    samp = get_cellassign_sample(adata, all_sigs[y])
    adata_filts[y] = samp
    save_file = 'CellAssignBI/nonLoggedCellAssignFiltST'+y+".h5ad"
    samp.write(save_file)


# In[ ]:


for y in adata_filts:    
    sample = adata_filts[y]
    scvi.external.CellAssign.setup_anndata(sample, size_factor_key="size_factor")
    gene_matrix =  create_marker_gene_matrix(all_sigs[y], list(sample.var.index))

    model = scvi.external.CellAssign(sample,gene_matrix)
    model.train()
    predictions = model.predict()
    write_file = 'CellAssignBI/CellAssignOutputST' +y+".h5ad"
    write_file_csv = 'CellAssignBI/CellAssignOutputST'+y+".csv"

    predictions.to_csv(write_file_csv)
    sample.write(write_file)
    
    
all_sigs = {
    "short": covid_sigs,
    "dev": laughney_sigs,
    "long": epi_sigs
}

def create_marker_gene_matrix(sigs_dict , sample_genes = []):
    if (len(sample_genes) == 0):
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
    else:
        all_genes = []
        for x in sigs_dict:
            all_genes.extend(sigs_dict[x])
        all_genes = list(set(all_genes) & set(sample_genes))
        df = []
        for x in sigs_dict:
            df_col = []
            for y in all_genes:
                df_col.append(int( y in sigs_dict[x]))
            df.append(df_col)
        return pd.DataFrame(df, index =list(sigs_dict.keys()), columns = all_genes).T

create_marker_gene_matrix(covid_sigs)


def get_cellassign_sample(sample, matrix):
    gene_matrix =  create_marker_gene_matrix(matrix)
    sample= sample[:,sample.var.index.isin(set(gene_matrix.index))]
    sample = sample.copy()
    return sample


# In[3]:


adata = sc.read_h5ad("filtered_cellassign_adata_BI.h5ad")



adata_filts = {}
for y in list(all_sigs.keys()):
    samp = get_cellassign_sample(adata, all_sigs[y])
    adata_filts[y] = samp
    save_file = 'CellAssignBI/nonLoggedCellAssignFilt'+y+".h5ad"
    samp.write(save_file)


# In[ ]:


for y in adata_filts:    
    sample = adata_filts[y]
    scvi.external.CellAssign.setup_anndata(sample, size_factor_key="size_factor")
    gene_matrix =  create_marker_gene_matrix(all_sigs[y], list(sample.var.index))

    model = scvi.external.CellAssign(sample,gene_matrix)
    model.train()
    predictions = model.predict()
    write_file = 'CellAssignBI/CellAssignOutput' +y+".h5ad"
    write_file_csv = 'CellAssignBI/CellAssignOutput'+y+".csv"

    predictions.to_csv(write_file_csv)
    sample.write(write_file)





