# Description: Integration and benchmarking for CRC cohort 
# Info: ~260k cells; 5000 HVGs; 59 batches | 1.5T mem; 120 cores
# Note: 3 integrate methods (Harmony, scVI, Seurat) with customized parameter benchmarking

# pkg
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import os
import re
import harmonypy as hm
import scanpy.external as sce
import itertools

sc.settings.verbosity = 3

os.makedirs('./data', exist_ok=True)
os.makedirs('./fig', exist_ok=True)
os.makedirs('./data/integrated', exist_ok=True)

# load data
adata = sc.read_h5ad('../01_QC/data/merged_CB_Filtered_data.h5ad.gz')
adata
# save the counts layer
adata.layers["counts"] = adata.X.copy()
# normalize to depth 10 000
sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
# logaritmize
sc.pp.log1p(adata)
# save normalized counts
adata.raw = adata
# gene selection
sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="cell_ranger", batch_key="SampleID", n_bins=20)
# check shared HVG
n_batches = adata.var["highly_variable_nbatches"].value_counts()
ax = n_batches.plot(kind="bar")
n_batches
# subset to 5k hvg
adata_hvg_5k = adata[:, adata.var["highly_variable"]].copy()
adata_hvg_5k
# save
save_file = './data/5k_hvg_unintegrate_data.h5ad.gz'
adata_hvg_5k.write_h5ad(save_file, compression='gzip')


# unintegrate obj
adata_before = adata_hvg_5k.copy()
# scale
sc.pp.scale(adata_before, max_value=10)
# pca
sc.tl.pca(adata_before, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_before, log=True, n_pcs=50)
# find nbhd
sc.pp.neighbors(adata_before n_neighbors=25, n_pcs=50)
sc.tl.umap(adata_before, min_dist=0.25, spread=1)
### plot
## sample 
# do things like this:
with plt.rc_context({"figure.figsize": (7, 7), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    sc.pl.umap(adata_before, color='SampleID', title='', size=10, alpha=0.5, show=False)
# or
with plt.rc_context({"figure.figsize": (7, 7), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    ax = sc.pl.umap(adata_before, size=5, show=False)
    sc.pl.umap(adata_before[adata_before.obs['Treatment'] == 'post_SCRT'], color='SampleID', title='', size=5, alpha=0.5, ax=ax, show=False)
## marker
# large cluster marker dict
marker_genes = {
    "T": ["CD3D", "CD3E", "CD3G"],
    "B": ["CD79A", "MS4A1", "CD79B"],
    "Plasma": ["SDC1", "MZB1", "IGHA2", "IGHG1", "IGHG3"],
    "Mast": ["KIT", "MS4A2", "TPSB2", "CPA3"],
    "Myeloid": ["CD1C", "CD1E", "CD68", "APOE", "CD163", "FCN1", "FCGR3B", "CD14"],
    "Stromal": ["COL1A1", "COL1A2", "ACTA2", "MYH11", "CLDN5", "VWF", "COL3A1"],
    "Epithelial": ["EPCAM", "KRT18", "KRT19", "CDH1"]
}

sc.pl.umap(adata_before, color=marker_genes['Epithelial'])



# Harmony integration
# PCA for harmony
adata_harmony = adata_hvg_5k.copy()
sc.pp.scale(adata_harmony, max_value=10)
sc.tl.pca(adata_harmony, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_harmony, log=True, n_pcs=50)
# save the unintegrated harmony obj (to avoid kernal carsh, `sce.pp.harmony_integrate` may cause jupyter crash..)
save_file = './data/adata_harmony_unintegrate_without_C007_BL.h5ad.gz'
adata_harmony.write_h5ad(save_file, compression='gzip')


############### 
####test1 ####: make harmony iter to convergence
###############
# Default max iter is 10
sce.pp.harmony_integrate(adata_harmony, key='SampleID', basis='X_pca',
                         adjusted_basis='X_pca_harmony', max_iter_harmony=50, max_iter_kmeans=50)
# save harmony-PCA embed
harmony_pca = adata_harmony.obsm['X_pca_harmony']
np.save("./data/5k_hvg_X_pca_harmony_convergence.npy", harmony_pca)

# umap
sc.pp.neighbors(adata_harmony, n_neighbors=25, n_pcs=50, use_rep='X_pca_harmony', key_added='harmony_neighbors')
sc.tl.umap(adata_harmony, min_dist=0.25, spread=1, neighbors_key='harmony_neighbors')
# save harmony_integrated anndata
save_file = './data/integrated/harmony_convergence/adata_harmony_convergence_integrated.h5ad.gz'
adata_harmony.write_h5ad(save_file, compression='gzip')

############### 
####test2 ####: make harmony do 50 iters
###############
sce.pp.harmony_integrate(adata_harmony, key='SampleID', basis='X_pca',
                         adjusted_basis='X_pca_harmony',
                         max_iter_harmony=50, max_iter_kmeans=50,
                         epsilon_cluster=float('-inf'), epsilon_harmony=float('-inf')
                         )
# save harmony-PCA embed
harmony_pca = adata_harmony.obsm['X_pca_harmony']
np.save("./data/5k_hvg_X_pca_harmony_50_iters.npy", harmony_pca)
# umap
sc.pp.neighbors(adata_harmony, n_neighbors=25, n_pcs=50, use_rep='X_pca_harmony', key_added='harmony_neighbors')
sc.tl.umap(adata_harmony, min_dist=0.25, spread=1, neighbors_key='harmony_neighbors')
# save harmony_integrated anndata
save_file = './data/integrated/harmony_50_iter/adata_harmony_50_iter_integrated.h5ad.gz'
adata_harmony.write_h5ad(save_file, compression='gzip')

# stop at here, reload `adata_hvg_5k`
adata_hvg_5k = sc.read_h5ad('./data/5k_hvg_unintegrate_data.h5ad.gz')


# scVI integration
# pkg
import scanpy as sc
import scvi
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scvi
import ray
import hyperopt
from ray import tune
from scvi import autotune

adata_scVI = adata_hvg_5k.copy()

############### 
####test1 ####: default
###############
adata_scVI = adata_hvg_5k.copy()
# default
scvi.model.SCVI.setup_anndata(adata_scVI,
                              layer="counts",
                              continuous_covariate_keys=["pct_counts_mt"],
                              batch_key="SampleID")
model_default = scvi.model.SCVI(adata_scVI)
model_default.train()
model_default.save("./test_scvi/scVI_5k_HVG_default/")

# load the model and get latent
model_default = scvi.model.SCVI.load("./test_scvi/scVI_5k_HVG_default/", adata=adata_scVI)
latent_default = model_default.get_latent_representation()
adata_scVI.obsm["X_scVI_default"] = latent_default
# umap
sc.pp.neighbors(adata_scVI, use_rep="X_scVI_default", n_neighbors=25, n_pcs=None, key_added='scVI_neighbors')
sc.tl.umap(adata_scVI, min_dist=0.25, spread=1, neighbors_key='scVI_neighbors')
# save scVI_default anndata
save_file = './integrated/scVI_default/adata_scVI_default_integrated.h5ad.gz'
adata_scVI.write_h5ad(save_file, compression='gzip')

############### 
####test2 ####: empirical argument
###############
adata_scVI = adata_hvg_5k.copy()
scvi.model.SCVI.setup_anndata(adata_scVI,
                              layer="counts",
                              continuous_covariate_keys=["pct_counts_mt"],
                              batch_key="SampleID")
# author empirical arguments
model_author = scvi.model.SCVI(adata_scVI, n_layers=2, n_latent=30, gene_likelihood="nb")
# check model
model_author.view_anndata_setup()
model_author
# training and save
model_author.train(max_epochs=300, early_stopping=True)
model_author.save("./test_scvi/scVI_5k_HVG_2/")
model_author = scvi.model.SCVI.load("./test_scvi/scVI_5k_HVG_2/", adata=adata_scVI)
# check loss
plt.plot(model_author.history['elbo_train']['elbo_train'], label='ELBO_train')
plt.plotmodel_author.history['elbo_validation']['elbo_validation'], label='ELBO_validation')
plt.legend()
plt.plot(model.history['reconstruction_loss_train']['reconstruction_loss_train'], label='train')
plt.plot(model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label='validation')
plt.legend()
# get letent and ump
latent_author = model_author.get_latent_representation()
adata_scVI.obsm["X_scVI_author"] = latent_author
sc.pp.neighbors(adata_scVI, use_rep="X_scVI_author", n_neighbors=25, n_pcs=None, key_added='scVI_neighbors')
sc.tl.umap(adata_scVI, min_dist=0.25, spread=1, neighbors_key='scVI_neighbors')
# save
save_file = './integrated/scVI_author/adata_scVI_author_integrated.h5ad.gz'
adata_scVI.write_h5ad(save_file, compression='gzip')

############### 
####test3 ####: ray tune
###############
# use ray tune
adata_hvg_5k = sc.read_h5ad('./data/5k_hvg_unintegrate_data.h5ad.gz')
adata_scVI = adata_hvg_5k.copy()
# setup
model_cls = scvi.model.SCVI
model_cls.setup_anndata(adata_scVI, layer="counts", continuous_covariate_keys=["pct_counts_mt"])
# define a new ModelTuner
scvi_tuner = autotune.ModelTuner(model_cls)
scvi_tuner.info()
# sarch space
search_space = {
    "n_latent": tune.choice([25, 30, 40, 50, 60]),
    "n_hidden": tune.choice([128, 192, 256]),
    "n_layers": tune.choice([1, 2, 3]),
    "lr": tune.loguniform(1e-4, 1e-2),
    "gene_likelihood": tune.choice(["nb", "zinb"])
}
ray.init(log_to_driver=False)
# do tune
results = scvi_tuner.fit(
    adata_scVI,
    metric="validation_loss",
    search_space=search_space,
    searcher="hyperopt",
    scheduler="asha",
    num_samples=40,
    max_epochs=150,
    resources={"cpu": 24, "gpu": 4},
)
# get best reuslt
results
print(results.model_kwargs)
print(results.train_kwargs)
# load all result
experiment_path = "./ray/tune_scvi_2023-05-03-01:01:39/"
print(f"Loading results from {experiment_path}...")
restored_tuner = tune.Tuner.restore(experiment_path)
result_grid = restored_tuner.get_results()
# check error
if result_grid.errors:
    print("One of the trials failed!")
else:
    print("No errors!")
# all result dataframe
result_grid.get_dataframe()
# plot
sns.boxplot(data=result_grid.get_dataframe(), y="validation_loss", x="config/n_latent", showfliers=False)
# plot all loss for each trails
with plt.rc_context({"figure.figsize": (5, 5), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    ax = None
    for result in result_grid:
        if result.config['n_latent'] == 60:
            label = f"{result.config['gene_likelihood']}; {result.config['n_hidden']}; {result.config['n_latent']}; {result.config['n_layers']}; {result.config['lr']:.5f}"
            if ax is None:
                ax = result.metrics_dataframe.plot("training_iteration", "validation_loss", label=label)
            else:
                result.metrics_dataframe.plot("training_iteration", "validation_loss", ax=ax, label=label)
    ax.set_title("Validation Loss vs. Training Iteration for All Trials")
    ax.set_ylabel("Validation Loss")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)

# retraining scVI model
adata_scVI = adata_hvg_5k.copy()
scvi.model.SCVI.setup_anndata(adata_scVI,
                              layer="counts",
                              continuous_covariate_keys=["pct_counts_mt"],
                              batch_key="SampleID")
model_optimum = scvi.model.SCVI(adata_scVI,
                        n_layers=2, n_hidden=256, n_latent=40,
                        dispersion='gene', gene_likelihood='zinb')
early_stopping_kwargs = {"early_stopping_monitor": "elbo_validation",
                         "early_stopping_patience": 15,
                         "early_stopping_min_delta": 0.001}
model_optimum.train(max_epochs=300, early_stopping=True,
            plan_kwargs={"lr": 0.0025},
            **early_stopping_kwargs
           )
# check ELBO
plt.plot(model_optimum.history['elbo_train']['elbo_train'], label='ELBO_train')
plt.plot(model_optimum.history['elbo_validation']['elbo_validation'], label='ELBO_validation')
plt.legend()
# save model
model_optimum.save("./test_scvi/scVI_5k_HVG_optimum/")
# load the before model
model_optimum = scvi.model.SCVI.load("./test_scvi/scVI_5k_HVG_optimum/", adata=adata_scVI)
adata_scVI.obsm["X_scVI_tune"] = latent_optimum
# umap
sc.pp.neighbors(adata_scVI, use_rep="X_scVI_tune", n_neighbors=25, n_pcs=None, key_added='scVI_neighbors')
sc.tl.umap(adata_scVI, min_dist=0.25, spread=1, neighbors_key='scVI_neighbors')
# save
save_file = './integrated/scVI_ray_tune_optimum/adata_scVI_ray_tune_integrated.h5ad.gz'
adata_scVI.write_h5ad(save_file, compression='gzip')

# stop at here, reload `adata_hvg_5k`
adata_hvg_5k = sc.read_h5ad('./data/5k_hvg_unintegrate_data.h5ad.gz')


# Seurat integration
# pkg
import numpy as np
import pandas as pd
import scanpy as sc
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import rpy2.rinterface_lib.callbacks
pandas2ri.activate()
anndata2ri.activate()
# rpy2 magics
%load_ext rpy2.ipython

adata_seurat = adata_hvg_5k.copy()
adata_seurat.obs["batch"] = adata_seurat.obs["batch"].astype(str)
del adata_seurat.uns

######## code for jupyter cell ########
%%R -i adata_seurat
adata_seurat
%%R -i adata_seurat
seurat_obj <- as.Seurat(adata_seurat, counts = "counts", data = "X")
seurat_obj
# save seurat obj
saveRDS(seurat_obj, "./data/seurat_obj_unintegrate.rds")
%%R
batch_list <- SplitObject(seurat_obj, split.by = 'batch')
batch_list
# save rds to avoid integrate error
saveRDS(batch_list, "./data/seurat_batch_list.rds")
######## code for jupyter cell ########

############### 
####test1 ####: ref-based rPCA
###############

########## R code ###################### R code ###################### R code ############
# workDir and libPath
setwd("/chenfeilab/Gaux/onGoingProject/wjw/CRC_project/cellbenderFilter/02_SampleIntegration")
.libPaths()
# lib
library(Seurat)
library(future)
# threads and mem
plan("multisession", workers = 48)
# 800g
options(future.globals.maxSize = 800000 * 1024^2)
# load data
batch_list <- readRDS("./data/seurat_batch_list.rds")
seurat_obj <- readRDS("./data/seurat_obj_unintegrate.rds")
features <- rownames(seurat_obj)
# before rpca, you must do pca
batch_list <- lapply(X = batch_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# default to find all mutual space, it is prohibitively computationally expensive, 1.5T mem can not run successfully
# so use ref-based
anchors_rpca <- FindIntegrationAnchors(batch_list, anchor.features=features, 
                                       reduction="rpca", dims = 1:50, 
                                       k.anchor = 5,
                                       reference = c(21,27,31,38,44,47,54,58),
                                       verbose = TRUE)
# save anchor
saveRDS(anchors_rpca, file = "./data/seurat_rpca_dim_50_ref_based_anchor.rds")
# integrate
integrated_rpca <- IntegrateData(anchors_rpca, dims = 1:50)
saveRDS(integrated_rpca, file = "./data/seurat_rpca_dim_50_ref_based_integrated.rds")

############### 
####test2 ####: ref-based CCA
###############
rm(list=ls())
# load
batch_list <- readRDS("./data/seurat_batch_list.rds")
seurat_obj <- readRDS("./data/seurat_obj_unintegrate.rds")
features <- rownames(seurat_obj)
# do cca ref-based find anchor
anchors_cca <- FindIntegrationAnchors(batch_list, anchor.features=features, 
                                      reduction="cca", dims = 1:50, 
                                      k.anchor = 5,
                                      reference = c(21,27,31,38,44,47,54,58),
                                      verbose = TRUE)
# save anchor
saveRDS(anchors_cca, file = "seurat_rpca_dim_50_ref_based_anchor.rds")
# do integrate
integrated_cca <- IntegrateData(anchors_cca, dims = 1:50)
saveRDS(integrated_cca, file = "./data/seurat_cca_dim_50_ref_based_integrated.rds")
########## R code ###################### R code ###################### R code ############

# switch to python to change seurat integrated counts matrixs to anndata

############### 
####test1 ####: ref-based rPCA
###############
adata_seurat_rPCA = adata_seurat.copy()
######## code for jupyter cell ########
%%R
# read the rPCA integrated data
integrated_rPCA <- readRDS("./data/seurat_rpca_dim_50_ref_based_integrated.rds")
integrated_rPCA
%%R -o integrated_rPCA_expr
integrated_rPCA_expr <- GetAssayData(integrated_rPCA)
integrated_rPCA_expr <- integrated_rPCA_expr[rownames(seurat_obj), colnames(seurat_obj)]
integrated_rPCA_expr <- t(integrated_rPCA_expr)
######## code for jupyter cell ########

adata_seurat_rPCA.X = integrated_rPCA_expr
adata_seurat_rPCA.layers["seurat_rPCA"] = integrated_rPCA_expr
print(adata_seurat_rPCA)
adata_seurat_rPCA.X

# do pca knn and ump
sc.tl.pca(adata_seurat_rPCA, svd_solver='arpack')
sc.pp.neighbors(adata_seurat_rPCA, n_neighbors=25, n_pcs=50)
sc.tl.umap(adata_seurat_rPCA, min_dist=0.25, spread=1)
# save
save_file = './data/integrated/seurat_rPCA/adata_seurat_rPCA_integrated.h5ad.gz'
adata_seurat_rPCA.write_h5ad(save_file, compression='gzip')


############### 
####test2 ####: ref-based CCA
###############
adata_seurat_CCA = adata_seurat.copy()

######## code for jupyter cell ########
%%R
integrated_CCA <- readRDS("./data/seurat_cca_dim_50_ref_based_integrated.rds")
integrated_CCA
%%R -o integrated_CCA_expr
integrated_CCA_expr <- GetAssayData(integrated_CCA)
integrated_CCA_expr <- integrated_CCA_expr[rownames(seurat_obj), colnames(seurat_obj)]
integrated_CCA_expr <- t(integrated_CCA_expr)
######## code for jupyter cell ########

# get np.float64 format, about 22G too large data, it will lead too slow to operate this object
adata_seurat_CCA.X = integrated_CCA_expr
adata_seurat_CCA.layers["seurat_CCA"] = integrated_CCA_expr
print(adata_seurat_CCA)
adata_seurat_CCA.X

# do pca knn and ump
sc.tl.pca(adata_seurat_CCA, svd_solver='arpack')
sc.pp.neighbors(adata_seurat_CCA, n_neighbors=25, n_pcs=50)
sc.tl.umap(adata_seurat_CCA, min_dist=0.25, spread=1)
# save
save_file = './data/integrated/seurat_CCA/adata_seurat_CCA_integrated.h5ad.gz'
adata_seurat_CCA.write_h5ad(save_file, compression='gzip')



# Benchmarking for integration
# pkg
import scib
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from adjustText import adjust_text
# raw data
adata = sc.read_h5ad('../01_QC/data/merged_CB_Filtered_data.h5ad.gz')
adata = adata[adata.obs['SampleID'] != 'C007_BL'].copy()
sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True)
sc.pp.log1p(adata)

# load all data
adata_hvg_5k = sc.read_h5ad('./data/5k_hvg_unintegrate_data.h5ad.gz')
adata_harmony_50_iter = sc.read_h5ad('./data/integrated/harmony_50_iter/adata_harmony_50_iter_integrated.h5ad.gz')
adata_harmony_convergence = sc.read_h5ad('./data/integrated/harmony_convergence/adata_harmony_convergence_integrated.h5ad.gz')
adata_scVI_tune = sc.read_h5ad('./data/integrated/scVI_ray_tune_optimum/adata_scVI_ray_tune_integrated.h5ad.gz')
adata_scVI_default = sc.read_h5ad('./data/integrated/scVI_default/adata_scVI_default_integrated.h5ad.gz')
adata_scVI_author = sc.read_h5ad('./data/integrated/scVI_author/adata_scVI_author_integrated.h5ad.gz')
adata_seurat_rPCA = sc.read_h5ad('./data/integrated/seurat_rPCA/adata_seurat_rPCA_integrated.h5ad.gz')
adata_seurat_CCA = sc.read_h5ad('./data/integrated/seurat_CCA/adata_seurat_CCA_integrated.h5ad.gz')

# get the cell-label free metrics
# for the feture space
adata_hvg_5k_metrics = scib.metrics.metrics(adata, adata_hvg_5k,
                                            batch_key='batch', label_key='Treatment',
                                            type_='full', organism='human', n_cores=48,
                                            ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
adata_seurat_rPCA_metrics = scib.metrics.metrics(adata, adata_seurat_rPCA,
                                                 batch_key='batch', label_key='Treatment',
                                                 type_='full', organism='human', n_cores=48,
                                                 ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
adata_seurat_CCA_metrics = scib.metrics.metrics(adata, adata_seurat_CCA,
                                                batch_key='batch', label_key='Treatment',
                                                type_='full', organism='human', n_cores=48,
                                                ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
# embed layer for scib
adata_harmony_50_iter.obsm['X_emb'] = adata_harmony_50_iter.obsm['X_pca_harmony']
adata_harmony_50_iter_metrics = scib.metrics.metrics(adata, adata_harmony_50_iter,
                                                     batch_key='batch', label_key='Treatment',
                                                     embed='X_pca_harmony', organism='human',
                                                     type_="embed", n_cores=48,
                                                     ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
adata_harmony_convergence.obsm['X_emb'] = adata_harmony_convergence.obsm['X_pca_harmony']
adata_harmony_convergence_metrics = scib.metrics.metrics(adata, adata_harmony_convergence,
                                                         batch_key='batch', label_key='Treatment',
                                                         embed='X_pca_harmony', organism='human',
                                                         type_="embed", n_cores=48,
                                                         ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
adata_scVI_author.obsm['X_emb'] = adata_scVI_author.obsm['X_scVI_author']                                                         
adata_scVI_author_metrics = scib.metrics.metrics(adata, adata_scVI_author,
                                                 batch_key='batch', label_key='Treatment',
                                                 embed='X_scVI_author', organism='human',
                                                 type_="embed", n_cores=48,
                                                 ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
adata_scVI_default.obsm['X_emb'] = adata_scVI_default.obsm['X_scVI_default']
adata_scVI_default_metrics = scib.metrics.metrics(adata, adata_scVI_default,
                                                  batch_key='batch', label_key='Treatment',
                                                  embed='X_scVI_default', organism='human',
                                                  type_="embed", n_cores=48,
                                                  ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)
adata_scVI_tune.obsm['X_emb'] = adata_scVI_tune.obsm['X_scVI_tune']
adata_scVI_tune_metrics = scib.metrics.metrics(adata, adata_scVI_tune,
                                               batch_key='batch', label_key='Treatment',
                                               embed='X_scVI_tune', organism='human',
                                               type_="embed", n_cores=48,
                                               ilisi_=True, pcr_=True, cell_cycle_=True, hvg_score_=True)

# summary
# Concatenate metrics results
metrics = pd.concat(
    [adata_hvg_5k_metrics,
     adata_harmony_convergence_metrics, adata_harmony_50_iter_metrics,
     adata_scVI_default_metrics, adata_scVI_author_metrics, adata_scVI_tune_metrics],
    axis="columns"
)
metrics = metrics.set_axis(
    ["Unintegrated_5k_hvg", "Harmony_convergence", "Harmony_50_iters", "scVI_default", "scVI_author", "scVI_tune"],
    axis="columns"
)
metrics.to_csv('./data/integrated/metrics_p1.csv')

metrics = pd.concat(
    [adata_seurat_rPCA_metrics, adata_seurat_CCA_metrics],
    axis="columns"
)
metrics = metrics.set_axis(
    ["Seurat_rPCA", "Seurat_CCA"],
    axis="columns"
)
metrics.to_csv('./data/integrated/metrics_p2.csv')

metrics_p1 = pd.read_csv('./data/integrated/metrics_p1.csv', index_col=0)
metrics_p2 = pd.read_csv('./data/integrated/metrics_p2.csv', index_col=0)
all_metrics = pd.concat([metrics_p1, metrics_p2], axis="columns")
all_metrics = all_metrics.iloc[[4, 5, 10, 12], :]
# deal the Harmony scale caused metrics different
all_metrics.iloc[3, [1,2]] = 1
all_metrics = all_metrics.T
all_metrics = all_metrics.iloc[:, [0, 2, 1, 3]]

# visualization
all_metrics.style.background_gradient(cmap="RdPu")

# rm non-uniform metric and scale the metrics
metrics = all_metrics.drop(columns=["hvg_overlap"])
metrics_scaled = (metrics - metrics.min()) / (metrics.max() - metrics.min())
metrics_scaled["Batch"] = metrics_scaled[["PCR_batch", "iLISI"]].mean(axis=1)
metrics_scaled["Bio"] = metrics_scaled[["cell_cycle_conservation"]].mean(axis=1)
metrics_scaled.style.background_gradient(cmap="Blues")

# plot
metrics_scaled['Type'] = metrics_scaled.index
with plt.rc_context({"figure.figsize": (12, 12), 'font.family': 'Arial', 'axes.facecolor': 'white'}):
    fig, ax = plt.subplots()
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    sns.scatterplot(data=metrics_scaled, x='Batch', y='Bio', hue='Type',
                    palette=['#394867', '#E5BEEC', '#917FB3', '#D6E8DB', '#C1D0B5', '#99A98F', '#E06469', '#F2B6A0'], ax=ax)
    annlist = []
    for k, v in metrics_scaled[["Batch", "Bio"]].iterrows():
        ann = ax.annotate(
            k,
            v
        )
        annlist.append(ann)
    plt.axhline(y=1, xmin=-0.1, xmax=1.1, color='grey', linestyle='-.', alpha=0.6)
    plt.axvline(x=1, ymin=-0.1, ymax=1.1, color='grey', linestyle='-.', alpha=0.6)
    ax.set_xlabel('Batch(PC regression score + iLISI)')
    ax.set_ylabel('Bio(cell_cycle)')
    plt.legend([], [], frameon=False)
    adjust_text(annlist)



# umap 
sc.pp.neighbors(adata_harmony, n_neighbors=25, n_pcs=50, use_rep='X_pca_harmony', key_added='harmony_neighbors')
sc.tl.umap(adata_harmony, min_dist=0.25, spread=1, neighbors_key='harmony_neighbors')
adata_harmony.obsm['X_umap_pc50_d25_s1'] = adata_harmony.obsm['X_umap'].copy()

args_dict = {
    "n_nbhds": [20, 25, 30],
    "n_pcs": [40, 50],
    "min_dist": [0.1, 0.15, 0.2, 0.25],
    "spread": [1, 1.5, 2]
}

# combination of `n_neighbors` and `n_pcs`:
for nbhd in args_dict['n_nbhds']:
    for pc in args_dict['n_pcs']:
        print(f"Computing: harmony_nbhd{nbhd}_pc{pc}")
        sc.pp.neighbors(adata_harmony, n_neighbors=nbhd, n_pcs=pc, use_rep='X_pca_harmony',
                        key_added=f"harmony_nbhd{nbhd}_pc{pc}")

# combination of `min_dist` and `spread` to get UMAP
keys, values = zip(*args_dict.items())
for nbhd, pc, dist, spread in itertools.product(*values):
    nbhd_pc = f"harmony_nbhd{nbhd}_pc{pc}"
    umap = f"X_umap_nbhd{nbhd}_pc{pc}_d{dist}_s{spread}".replace('.', 'p')
    print("Computing " + umap + " for " + nbhd_pc)
    sc.tl.umap(adata_harmony, min_dist=dist, spread=spread, neighbors_key=nbhd_pc)
    print("Saving " + umap)
    adata_harmony.obsm[umap] = adata_harmony.obsm['X_umap'].copy()

# save all combine adata_obj
save_file = './data/adata_harmony_all_combine_without_C007_BL.h5ad.gz'
adata_harmony.write_h5ad(save_file, compression='gzip')

# plot
umapRegex = re.compile('X_umap_nbhd[0-9]{2}_pc[0-9]{2}_d0p[0-9]*_s[0-9]*p?[0-9]?')
all_umap = [umap for umap in list(adata_harmony.obsm.keys()) if umapRegex.search(umap)]

def return_plot(umap: str, axe):
    return sc.pl.embedding(
        adata_harmony,
        basis=umap,
        title=umap,
        color='SampleID',
        size=5,
        alpha=0.5,
        legend_loc=None,
        show=False,
        ax=axe
    )

plot_all_umap = all_umap
n_col = 6
n_row = 12

fig, ax_array = plt.subplots(n_row, n_col, squeeze=False, figsize=(n_col * 4, n_row * 4))
for i, ax_row in enumerate(ax_array):
    for j, (umap, axes) in enumerate(zip(plot_all_umap, ax_row)):
            # axes.set_title('{},{}'.format(i, j))
            return_plot(umap=umap, axe=axes)
            axes.set_xlabel('UMAP_1')
            axes.set_ylabel('UMAP_2')
            print('Ploting: ({},{}): {}'.format(i, j, umap))
            # new line to plot
            if (j+1) % n_col == 0:
                del all_umap[0:n_col]
plt.savefig("./fig/harmony_all_umap_combination.png", bbox_inches='tight')
plt.close()


