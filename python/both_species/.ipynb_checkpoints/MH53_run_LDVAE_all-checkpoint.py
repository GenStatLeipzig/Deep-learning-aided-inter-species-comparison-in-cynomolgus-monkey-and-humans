import numpy as np
import torch
import scanpy as sc
import matplotlib.pyplot as plt
#import torchmetrics
from torchmetrics.utilities.data import dim_zero_sum 
import sys
import os
path_helper = ["C:\\","Users","vfriedrich","projects","monkey_IZI","git_documentation","scRNAseq_cross_species_primate_human","analysis","helper"]
sys.path.append(os.path.join(*path_helper))
import helperVDF as h
print(sys.executable)
import scarches
import anndata
import pandas as pd
import anndata as ad
import scvi
import gzip
from sklearn.model_selection import train_test_split

pre = "MH53"
drive = 'F'
base_model_path,base_table_path,base_plots_path,base_anndata_objects = h.return_local_paths(drive = drive,
                                                                                            pre = pre,
                                                                                            add_path = True)

_,_,_,base_anndata_objects_MH44 = h.return_local_paths(drive = drive,pre = "MH44",add_path = True)

adata = sc.read_h5ad(os.path.join(base_anndata_objects_MH44,'MH44_clusterannoVAE.h5ad'))
adata.obs.rename(columns={"timepoint_x": "timepoint", "individual_x": "individual"}, inplace=True)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata_train = adata[adata.obs['individual'].isin(['cyno1','human1'])]
adata_test = adata[adata.obs['individual'].isin(['cyno2','human2'])]
adata_train.write(os.path.join(base_anndata_objects,pre + '_adata_train.h5ad'))
adata_test.write(os.path.join(base_anndata_objects,pre + '_adata_test.h5ad'))

celltypes = ['CD4 T','NK','CD8 T','B','CD14 Mono','CD16 Mono','MAIT','NK+Proliferating']
celltype_col = 'cluster_azimut1_5_scanvi'
n_latent = 12
batch_size = 32

for celltype in celltypes:
    if celltype != 'NK+Proliferating':
        adata_ct = h.filter_adata_obs(adata=adata_train,col_name=celltype_col,val=celltype)
    else:
        adata_NK = h.filter_adata_obs(adata=adata_train,col_name=celltype_col,val='NK')
        adata_NK_prolif = h.filter_adata_obs(adata=adata_train,col_name=celltype_col,val='NK Proliferating')
        adata_ct = ad.concat([adata_NK,adata_NK_prolif])
    adata_ct = adata_ct.copy()
    scvi.model.LinearSCVI.setup_anndata(adata_ct)
    model = scvi.model.LinearSCVI(adata_ct, n_latent=n_latent)
    model.train(max_epochs=100,train_size=0.9,batch_size=batch_size,early_stopping=True,early_stopping_min_delta=0.001,early_stopping_patience=25,plan_kwargs={"lr": 5e-3})
    save_string = os.path.join(base_model_path,pre + '_LDVAE_model_' + str(celltype) +'.pt')
    model.save(save_string,save_anndata=True)


base_package_version_path = h.return_package_version_local_path(drive=drive)
h.save_package_versions(base_package_version_path,pre,do_print = True)
h.print_main_versions()