from anndata import AnnData
import anndata
from scipy import sparse, io
import scipy
import pandas as pd
import scipy.io
import os
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
matplotlib.use('TkAgg')
import numpy as np
import seaborn as sb
import math
import scanpy.external as sce
#import scrublet as scr

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (7,7)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")

###############################################################################################
#################### Only Vascular Smooth muscle cells
###############################################################################################
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri
pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython
%%R
library(scran)
library(dplyr)






%config InlineBackend.figure_format = 'retina'

test3_fb = anndata.AnnData(X=test3[test3.obs['leiden_r05'].isin(['3', '6', '9'])].layers['counts'], obs=test3[test3.obs['leiden_r05'].isin(['3', '6', '9'])].obs, var=test3[test3.obs['leiden_r05'].isin(['3', '6', '9'])].var)
test3_fb.layers["counts"] = test3_fb.X.copy()

adata_pp = test3_fb.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_fb.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_fb.obs['size_factors'] = size_factors

test3_fb.X /= test3_fb.obs['size_factors'].values[:, None]
test3_fb.X = scipy.sparse.csr_matrix(test3_fb.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)

test3_fb.layers['scran'] = test3_fb.X

sc.pp.log1p(test3_fb) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_fb.layers['scran_log1p'] = test3_fb.X

test3_fb.raw = test3_fb ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

sc.pp.highly_variable_genes(test3_fb)

test3_fb.var['highly_variable'].value_counts() # 2,612 ==> 2021-08-20, # 2,941 ==> 2021-09-28

sc.pp.filter_genes(test3_fb, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE

sc.pp.scale(test3_fb, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
# adata.raw.X의 mean 과 std를 output함
sc.tl.pca(test3_fb, n_comps=100, use_highly_variable=True, svd_solver='arpack')
#sc.pl.pca_variance_ratio(test3_fb, n_pcs=100)

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_fb, batch_key='batch', n_pcs=10, neighbors_within_batch=5, trim=None) #####
sc.tl.umap(test3_fb, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
#test3_endo.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_fb, resolution=0.5, key_added='fb_leiden_r05')
sc.tl.leiden(test3_fb, resolution=1.0, key_added='fb_leiden_r10')

test3_fb.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
import magic
test3_vsmc_MAGIC = test3_fb.copy()
test3_vsmc_MAGIC.X = test3_fb.layers['scran_log1p']
test3_vsmc_MAGIC = magic.MAGIC().fit_transform(test3_vsmc_MAGIC)
test3_fb.layers['magic'] = test3_vsmc_MAGIC.X

sc.pl.umap(test3_fb, color=['fb_leiden_r05', 'fb_leiden_r10', 'leiden_r05', 'batch'], add_outline=False, legend_loc='right margin', size=50, color_map=cmap, ncols=4)

