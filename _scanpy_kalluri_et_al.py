#export PATH=/data/Projects/phenomata/99.Tools/anaconda3/bin:$PATH
#source activate scanpy_1.8.1

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
import seaborn as sns
import math
import scanpy.external as sce
#import scrublet as scr
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (6,6)
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Arial'
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
batch_palette=['#689aff', '#fdbf6f', '#b15928']
%matplotlib
%autoindent

#test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
#test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

m03_1 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Kalluri_et_al_2019/WT1_1.h5")
m03_1.var_names_make_unique()

m03_2 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Kalluri_et_al_2019/WT1_2.h5")
m03_2.var_names_make_unique()

m03_3 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Kalluri_et_al_2019/WT3_1.h5")
m03_3.var_names_make_unique()

m03_4 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Kalluri_et_al_2019/WT3_2.h5")
m03_4.var_names_make_unique()

mito_genes = m03_1.var_names.str.startswith('mt-')
m03_1.obs['percent_mito'] = np.ravel(np.sum(m03_1[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m03_1.X, axis=1))
m03_2.obs['percent_mito'] = np.ravel(np.sum(m03_2[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m03_2.X, axis=1))
m03_3.obs['percent_mito'] = np.ravel(np.sum(m03_3[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m03_3.X, axis=1))
m03_4.obs['percent_mito'] = np.ravel(np.sum(m03_4[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m03_4.X, axis=1))
sce.pp.scrublet(m03_1, adata_sim=None, sim_doublet_ratio=2.0, expected_doublet_rate=0.05, stdev_doublet_rate=0.02, synthetic_doublet_umi_subsampling=1.0, knn_dist_metric='euclidean', n_prin_comps=30, verbose=True)
sce.pp.scrublet(m03_2, adata_sim=None, sim_doublet_ratio=2.0, expected_doublet_rate=0.05, stdev_doublet_rate=0.02, synthetic_doublet_umi_subsampling=1.0, knn_dist_metric='euclidean', n_prin_comps=30, verbose=True)
sce.pp.scrublet(m03_3, adata_sim=None, sim_doublet_ratio=2.0, expected_doublet_rate=0.05, stdev_doublet_rate=0.02, synthetic_doublet_umi_subsampling=1.0, knn_dist_metric='euclidean', n_prin_comps=30, verbose=True)
sce.pp.scrublet(m03_4, adata_sim=None, sim_doublet_ratio=2.0, expected_doublet_rate=0.05, stdev_doublet_rate=0.02, synthetic_doublet_umi_subsampling=1.0, knn_dist_metric='euclidean', n_prin_comps=30, verbose=True)

#integrated = AnnData.concatenate(m03_1, m03_2, m03_3, m03_4, join='outer', batch_categories = ['m03_1', 'm03_2', 'm03_3', 'm03_4'], index_unique = '-')
#sc.pp.calculate_qc_metrics(integrated, log1p=False, inplace=True)
#sns.scatterplot(data=integrated.obs[['batch', 'n_genes_by_counts', 'total_counts', 'percent_mito']], x="total_counts", y="n_genes_by_counts", hue="percent_mito")
#sns.histplot(data=integrated.obs, x='percent_mito', hue='batch', kde=True, stat='density', common_bins=False, common_norm=False)
#sc.pl.scatter(integrated, x='total_counts', y='n_genes_by_counts', color='percent_mito')

sc.pp.filter_cells(m03_1, min_genes=2000)
sc.pp.filter_cells(m03_2, min_genes=2000)
sc.pp.filter_cells(m03_3, min_genes=2000)
sc.pp.filter_cells(m03_4, min_genes=2000)

m03_1 = m03_1[m03_1.obs['percent_mito'] < 0.2]
m03_2 = m03_2[m03_2.obs['percent_mito'] < 0.2]
m03_3 = m03_3[m03_3.obs['percent_mito'] < 0.2]
m03_4 = m03_4[m03_4.obs['percent_mito'] < 0.2]


integrated = AnnData.concatenate(m03_1, m03_2, m03_3, m03_4, join='outer', batch_categories = ['m03_1', 'm03_2', 'm03_3', 'm03_4'], index_unique = '-')

sc.pp.filter_genes(integrated, min_cells=5) # integrated.var??? n_cells ??????
integrated.layers["counts"] = integrated.X.copy()
integrated.raw = integrated

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

adata_pp = integrated.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp) # works on anndata.X
sc.tl.pca(adata_pp, n_comps=15) ## ????????? ??? n_component??? ????????? ????????? size_factors??? estimation?????? ??? ????????? ???????
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = integrated.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
integrated.obs['size_factors'] = size_factors

integrated.X /= integrated.obs['size_factors'].values[:, None]
integrated.layers['scran'] = integrated.X # For cellphoneDB
sc.pp.log1p(integrated) # works on anndata.X
integrated.X = scipy.sparse.csr_matrix(integrated.X)
integrated.raw = integrated ## ==> log transforamtion ??? ?????? raw??? ???????????? ???.

test3 = integrated.copy()
test3.raw = test3
test3.layers['scran_log1p'] = test3.X

sc.pp.highly_variable_genes(test3)
test3.var['highly_variable'].value_counts() # 2,410 ==> 2021-08-10 # 2,513 ==>

sc.pp.scale(test3, max_value=10) # tabula muris senis default (2021-08-10) # mean and std on adata.var
#sc.pp.scale(test3, zero_center=True, max_value=10, copy=False, layer=None, obsm=None)
cell_cycle_genes=[x.strip()[0] + x.strip()[1:].lower() for x in open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/regev_lab_cell_cycle_genes.txt")]
s_genes= cell_cycle_genes[:43]
g2m_genes= cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in test3.var_names]
sc.tl.score_genes_cell_cycle(test3, s_genes=s_genes, g2m_genes=g2m_genes)

sc.tl.pca(test3, n_comps=100, use_highly_variable=True, svd_solver='arpack')
#matplotlib.use('TkAgg')
#%matplotlib
#sc.pl.pca_variance_ratio(test3, n_pcs=50, log=True)
#sc.pl.pca(test3, color=['batch'], legend_loc='right margin', size=8, add_outline=False, color_map='CMRmap', components=['1,2'])

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3, batch_key='batch', n_pcs=15, neighbors_within_batch=20, trim=None) #####

#sc.pp.neighbors(test3, n_neighbors=15, n_pcs=15, method='umap')
sc.tl.umap(test3, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')

#test3.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
test3.obs['Doublet'] = test3.obs['predicted_doublet'].astype(str).astype('category')
sc.tl.leiden(test3, resolution=0.5, key_added='leiden_r05') #### 0 ~ 13 ==> 2021-09-28
sc.tl.leiden(test3, resolution=1.0, key_added='leiden_r10')
sc.pl.umap(test3, color=['batch', 'leiden_r05', 'leiden_r10', 'Cdh5'], add_outline=False, legend_loc='right margin', size=100, color_map=cmap)

sc.tl.rank_genes_groups(test3, 'leiden_r05', method='wilcoxon', corr_method='benjamini-hochberg', use_raw=True, pts=True) # key_added=''
sc.pl.rank_genes_groups(test3, n_genes=5, sharey=False)
sc.pl.rank_genes_groups_heatmap(test3, n_genes=10, min_logfoldchange=2, cmap='cividis', show_gene_labels=True)

markers = ["Pecam1", "Cdh5", "Nos3", "Acta2", "Cnn1", "Tagln", "Rgs5", "Kcnj8", "Col1a1", "Col5a1", "Dpt", "Cd19", "Ighm", "Cd14", "Cd68", "Cd3d"] # Cd3g ??????
sc.pl.stacked_violin(test3, markers, groupby='batch')

result = test3.uns['rank_genes_groups']
groups = result['names'].dtype.names
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})
deg_wilcoxon.to_csv("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/20210916_scanpy_deg.csv", mode='w')

# leiden_r0.5 ==> 20210916_scanpy_deg.csv ==> Fold change cutoff 1.0 / padj < 0.05 ==> canonical markers (2021-11-16)
leiden_to_celltype_dict = {'0': 'Smooth muscle cells',
'1': 'Smooth muscle cells',
'2': 'Smooth muscle cells',
'3': 'Fibroblasts',
'4': 'Smooth muscle cells',
'5': 'Endothelial cells',
'6': 'Fibroblasts',
'7': 'Endothelial cells',
'8': 'Smooth muscle cells',
'9': 'Fibroblasts',
'10': 'B cells',
'11': 'M\u03A6',
'12': 'T cells',
'13': 'Fibroblasts'}
test3.obs['celltype'] = test3.obs['leiden_r05'].map(lambda x: leiden_to_celltype_dict[x]).astype('category')

sc.pl.umap(test3, color=['batch', 'leiden_r05', 'celltype'], add_outline=False, legend_loc='right margin', size=50, color_map='CMRmap')

celltype_marker = {'EC': ['Pecam1', 'Cdh5', 'Vwf', 'Nos3'],
'SMC': ['Acta2', 'Tagln', 'Cnn1', 'Cnn2'],
'FB': ['Dpt', 'Col1a1', 'Col5a1', 'Pdgfra'],
'B cells': ['Ighm', 'Cd19'],
'M??':['Cd14', 'Cd68'],
'T cells':['Cd3d', 'Cd3g']}
reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'M\u03A6', 'T cells')
#reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'Macrophages', 'T cells')
test3.obs['celltype'] = test3.obs['celltype'].cat.reorder_categories(list(reordered), ordered=True)

sc.pl.umap(test3, color=['celltype'], add_outline=False, legend_loc='right margin', size=30, color_map=cmap, palette='Paired')

leiden_to_celltype2_dict = {'0': 'VSMC_1',
'1': 'VSMC_2',
'2': 'VSMC_3',
'3': 'FB_1',
'4': 'VSMC_4',
'5': 'EC_1',
'6': 'FB_2',
'7': 'EC_2',
'8': 'VSMC_5',
'9': 'FB_3',
'10': 'B-lympho',
'11': 'M\u03A6',
'12': 'T-lmpho',
'13': 'Neuronal'}
test3.obs['celltype2'] = test3.obs['leiden_r05'].map(lambda x: leiden_to_celltype2_dict[x]).astype('category')
sc.pl.umap(test3, color=['celltype2'], add_outline=False, legend_loc='on data', size=30, color_map=cmap, palette='tab20')

colormap = dict()
c = 0
for i in leiden_to_celltype2_dict.values():
    colormap[i] = test3.uns['celltype2_colors'][c]
    c += 1

df = pd.concat([test3.obs['batch'], test3.obs['celltype2']], axis=1)
ax = pd.crosstab(df['batch'], df['celltype2'], normalize=0).plot.bar(stacked=True, color=colormap)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()
