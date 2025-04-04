import torchmetrics
import scgen
import numpy as np
import scanpy as sc
import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.spatial.distance import cdist

def mean_euclidean_distance(embedding1,embedding2):
    return np.mean(np.linalg.norm(embedding1 - embedding2,axis=1))

import warnings

import sys
print(sys.executable)
cluster_pre = "MH49_c"
pre = "MH66_c"

#feature importance analysis for models from MH29_v2_c
base_model_path = '/work/username/output/MH49_c/models'
save_loc_tables = '/work/username/output/MH66_c/tables'

celltypes = ['B',
 'CD14 Mono',
 'CD16 Mono',
 'CD4 T',
 'CD8 T',
 'MAIT',
 'NK+Proliferating']

columns_lat_dim = [f'latent_dim_{i+1}' for i in range(10)]
columns_df =  columns_lat_dim + ['celltype','gene']
for celltype in celltypes:
    loader = os.path.join(base_model_path,celltype,cluster_pre + '_scgen_model_' + celltype +  '.pt')
    model = scgen.SCGEN.load(loader)
    model.adata.obs.rename(columns={'species_x': 'species'}, inplace=True)
    genes = list(model.adata.var.index)
    full_embedding = model.get_latent_representation(model.adata)
    df_perturb = pd.DataFrame(columns= columns_df)
    for gene in genes:
        adata = model.adata.copy()
        full_embedding = model.get_latent_representation(model.adata)
        index = adata.var.index.tolist().index(gene)
        
        human_cynosk = adata.obs['species'] == 'human'
        cyno_cynosk = adata.obs['species'] == 'cyno'
        mean_val_human = np.mean(adata.X[human_cynosk,index])
        var_val_human = np.var(adata.X[human_cynosk,index].toarray())
        mean_val_cyno = np.mean(adata.X[cyno_cynosk ,index])
        var_val_cyno = np.var(adata.X[cyno_cynosk ,index].toarray())
    
        adata.obs['species_mean'] = pd.Series(index=adata.obs.index, dtype=str)
        adata.obs['species_var'] = pd.Series(index=adata.obs.index, dtype=str)


        adata.obs.loc[adata.obs['species'] == 'human', 'species_mean'] = mean_val_human
        adata.obs.loc[adata.obs['species'] == 'cyno', 'species_mean'] = mean_val_cyno
        adata.obs.loc[adata.obs['species'] == 'human', 'species_var'] = var_val_human
        adata.obs.loc[adata.obs['species'] == 'cyno', 'species_var'] = var_val_cyno

        adata.obs['counts_input'] = adata.X[:,index].toarray().flatten().copy()
        adata.obs['diff_to_mean'] = (adata.obs['counts_input'] - adata.obs['species_mean']) 
        # if variance < 1 per species - shift to mean, otherwise shift negative proportional to variance towards mean
        adata.obs['species_var_pert_factor'] = np.where(adata.obs['species_var'] < 1, 1, adata.obs['species_var'])
        adata.obs['species_var_pert_factor'] = 1/adata.obs['species_var_pert_factor']
        adata.obs['perturb_vector'] = adata.obs['diff_to_mean']*adata.obs['species_var_pert_factor']
        adata.obs['perturbed_input'] = adata.obs['counts_input'] - adata.obs['perturb_vector']
        adata.X[:,index] = adata.obs['perturbed_input'].values
        lat_rep = model.get_latent_representation(adata)
        row_iter = pd.DataFrame([np.linalg.norm(lat_rep - full_embedding,axis=0)], columns=columns_lat_dim)
        row_iter['celltype'] = celltype
        row_iter['gene'] = gene
        df_perturb = df_perturb.append(row_iter,ignore_index=True)
    df_perturb.to_csv(os.path.join(save_loc_tables, pre + '_' + celltype + '_pert.csv'))
    
import pkg_resources
with open(os.path.join(os.getcwd(),'package_versions', pre + '_package_versions.txt'), "w") as file:
        for package in pkg_resources.working_set:
            if package.key != 'batchglm':
                file.write(f"{package.key}=={package.version}\n")
                print(f"{package.key}=={package.version}")
