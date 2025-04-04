import numpy as np
import torch
import scanpy as sc
import matplotlib.pyplot as plt
#import torchmetrics
from torchmetrics.utilities.data import dim_zero_sum 
import sys
import os
path_helper = ["C:\\","Users","username","projects","monkey_IZI","git_documentation","scRNAseq_cross_species_primate_human","analysis","helper"]
sys.path.append(os.path.join(*path_helper))
import helperVDF as h
print(sys.executable)
import scarches
import anndata
import pandas as pd
import scvi
import gzip
from sklearn.model_selection import train_test_split

pre = "MH46"
drive = 'F'
base_model_path,base_table_path,base_plots_path,base_anndata_objects = h.return_local_paths(drive = drive,
                                                                                            pre = pre,
                                                                                            add_path = True)

_,_,_,base_anndata_objects_MH44 = h.return_local_paths(drive = drive,pre = "MH44",add_path = True)

adata = sc.read_h5ad(os.path.join(base_anndata_objects_MH44,'MH44_clusterannoVAE.h5ad'))
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.subsample(adata,n_obs=1000,random_state=5)
train_indices, test_indices = train_test_split(adata.obs.index, test_size=0.2, random_state=5)

train_adata = adata[train_indices].copy() 
test_adata = adata[test_indices].copy()

train_adata.write(os.path.join(base_anndata_objects,pre + '_train_adata_hyperparameter.h5ad'))
test_adata.write(os.path.join(base_anndata_objects,pre + '_test_adata_hyperparameter.h5ad'))
scvi.model.LinearSCVI.setup_anndata(train_adata)

n_latents = [6,8,10,12]
batch_sizes = [8,16,32,64]
for n_latent in n_latents:
    for batch_size in batch_sizes:
        model = scvi.model.LinearSCVI(train_adata, n_latent=n_latent)
        model.train(max_epochs=100,train_size=0.9,batch_size=batch_size,early_stopping=True,early_stopping_min_delta=0.001,early_stopping_patience=25,plan_kwargs={"lr": 5e-3})
        save_string = os.path.join(base_model_path,pre + '_LDVAE_model_nlat' + str(n_latent) + '_bs' + str(batch_size)  + '.pt')
        model.save(save_string,save_anndata=True)

base_package_version_path = h.return_package_version_local_path(drive=drive)
h.save_package_versions(base_package_version_path,pre,do_print = True)
h.print_main_versions()