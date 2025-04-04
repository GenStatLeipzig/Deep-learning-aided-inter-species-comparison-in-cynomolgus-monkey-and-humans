import numpy as np
import scanpy as sc
import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import sys
#path_helper = '/home/username/compute/primates/MH29_c/helper'
#sys.path.append(os.path.join(*path_helper))
#import helperVDF as h
#import helper_scgen as hscg
#print(sys.executable)
import anndata as ad
#import time
import scgen
import scvi
from sklearn.model_selection import train_test_split

pre = 'MH45_c'


output_path = '/work/username/output'
input_path = '/work/username/input_data'
base_anndata_objects = '/work/username/output/'+ pre + '/anndata_objects'
base_model_path = '/work/username/output/'+ pre + '/models'
adata = sc.read_h5ad(os.path.join(input_path,'MH44_clusterannoVAE.h5ad'))

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.subsample(adata,n_obs=1000,random_state=5)
train_indices, test_indices = train_test_split(adata.obs.index, test_size=0.2, random_state=5)

train_adata = adata[train_indices].copy() 
test_adata = adata[test_indices].copy()

train_adata.write(os.path.join(base_anndata_objects,pre + '_train_adata_hyperparameter.h5ad'))
test_adata.write(os.path.join(base_anndata_objects,pre + '_test_adata_hyperparameter.h5ad'))
scvi.data.setup_anndata(train_adata)

n_latents = [6,8,10,12]
batch_sizes = [8,16,32,64]

for n_latent in n_latents:
    for batch_size in batch_sizes:
        model = scgen.SCGEN(train_adata,n_latent = n_latent)
        model.train(max_epochs=100,train_size=0.9,batch_size=batch_size,early_stopping=True,early_stopping_min_delta=0.001,early_stopping_patience=25)
        save_string = os.path.join(base_model_path,'nlat' +str(n_latent) + '_bs' +  str(batch_size),pre + '_scgen_model_nlat' + str(n_latent) + '_bs' + str(batch_size)  + '.pt')
        model.save(save_string,save_anndata=True)
