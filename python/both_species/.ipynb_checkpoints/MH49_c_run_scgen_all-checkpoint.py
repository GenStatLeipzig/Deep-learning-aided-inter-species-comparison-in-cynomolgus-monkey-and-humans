import numpy as np
import scanpy as sc
import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import sys
import anndata as ad
#import time
import scgen
import scvi
from sklearn.model_selection import train_test_split

def filter_adata_obs(adata,col_name,val):
    '''
    adata    :   anndata object
    col_name :   string, column name
    val      :   string or float, value in column to filter for
    '''
    return adata[adata.obs[col_name] == val]

pre = 'MH49_c'


output_path = '/work/username/output'
input_path = '/work/username/input_data'
base_anndata_objects = '/username/output/'+ pre + '/anndata_objects'
base_model_path = '/work/username/output/'+ pre + '/models'
adata = sc.read_h5ad(os.path.join(input_path,'MH44_clusterannoVAE.h5ad'))
adata.obs.rename(columns={"timepoint_x": "timepoint", "individual_x": "individual"}, inplace=True)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
celltypes = ['CD4 T','NK','CD8 T','B','CD14 Mono','CD16 Mono','MAIT','NK+Proliferating']
celltype_col = 'cluster_azimut1_5_scanvi'
n_latent = 10
batch_size = 8

for celltype in celltypes:
    if celltype != 'NK+Proliferating':
        adata_ct = filter_adata_obs(adata=adata,col_name=celltype_col,val=celltype)
    else:
        adata_NK = filter_adata_obs(adata=adata,col_name=celltype_col,val='NK')
        adata_NK_prolif = filter_adata_obs(adata=adata,col_name=celltype_col,val='NK Proliferating')
        adata_ct = ad.concat([adata_NK,adata_NK_prolif])
    adata_ct = adata_ct.copy()
    scvi.data.setup_anndata(adata_ct)
    model = scgen.SCGEN(adata_ct,n_latent = n_latent)
    model.train(max_epochs=100,train_size=0.9,batch_size=batch_size,early_stopping=True,early_stopping_min_delta=0.001,early_stopping_patience=25)
    print(model.history)
    save_string = os.path.join(base_model_path,celltype,pre + '_scgen_model_' + str(celltype) + '.pt')
    model.save(save_string,save_anndata=True)

