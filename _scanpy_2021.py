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
import math
import scanpy.external as sce
import scrublet as scr

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (7,7)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

m01 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/01.Cell-Ranger/01month_filtered_feature_bc_matrix.h5")
m01.var_names_make_unique()

m10 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/01.Cell-Ranger/10months_filtered_feature_bc_matrix.h5")
m10.var_names_make_unique()

m20 = sc.read_10x_h5("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/01.Cell-Ranger/20months_filtered_feature_bc_matrix.h5")
m20.var_names_make_unique()

mito_genes = m01.var_names.str.startswith('mt-')
m01.obs['percent_mito'] = np.ravel(np.sum(m01[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m01.X, axis=1))
m10.obs['percent_mito'] = np.ravel(np.sum(m10[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m10.X, axis=1))
m20.obs['percent_mito'] = np.ravel(np.sum(m20[:, mito_genes].X, axis=1)) / np.ravel(np.sum(m20.X, axis=1))

sc.pp.filter_cells(m01, min_counts=2000)
sc.pp.filter_cells(m01, min_genes=1500)

sc.pp.filter_cells(m10, min_counts=3000)
sc.pp.filter_cells(m10, min_genes=1500)

sc.pp.filter_cells(m20, min_counts=3000)
sc.pp.filter_cells(m20, min_genes=1500)

m01 = m01[m01.obs['percent_mito'] < 0.2]
m10 = m10[m10.obs['percent_mito'] < 0.2]
m20 = m20[m20.obs['percent_mito'] < 0.2]

integrated = AnnData.concatenate(m01, m10, m20, join='outer', batch_categories = ['m01', 'm10', 'm20'], index_unique = '-')

sc.pp.filter_genes(integrated, min_cells=5) # integrated.var에 n_cells 추가
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
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
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
integrated.raw = integrated ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

test3 = integrated.copy()
test3.raw = test3
test3.layers['scran_log1p'] = test3.X

sc.pp.highly_variable_genes(test3)
test3.var['highly_variable'].value_counts() # 2,410 ==> 2021-08-10 # 2,513 ==>

sc.pp.scale(test3, max_value=10) # tabula muris senis default (2021-08-10) # mean and std on adata.var

cell_cycle_genes=[x.strip()[0] + x.strip()[1:].lower() for x in open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/regev_lab_cell_cycle_genes.txt")]
s_genes= cell_cycle_genes[:43]
g2m_genes= cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in test3.var_names]


sc.tl.pca(test3, n_comps=100, use_highly_variable=True, svd_solver='arpack')
matplotlib.use('TkAgg')
%matplotlib
#sc.pl.pca_variance_ratio(test3, n_pcs=50, log=True)
#sc.pl.pca(test3, color=['batch'], legend_loc='right margin', size=8, add_outline=False, color_map='CMRmap', components=['1,2'])

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None) #####
sc.tl.umap(test3, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
#test3.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
test3.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
#sc.pl.umap(test3, color=['batch'], add_outline=False, legend_loc='right margin', size=20, color_map='CMRmap')

sc.tl.leiden(test3, resolution=0.5, key_added='leiden_r05') #### 0 ~ 13 ==> 2021-09-28
sc.tl.leiden(test3, resolution=1.0, key_added='leiden_r10')
sc.pl.umap(test3, color=['batch', 'leiden_r05', 'leiden_r10'], add_outline=False, legend_loc='right margin', size=20, color_map='CMRmap')

sc.tl.rank_genes_groups(test3, 'leiden_r05', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(test3, n_genes=5, sharey=False)
sc.pl.rank_genes_groups_heatmap(test3, n_genes=10, min_logfoldchange=2, cmap='cividis', show_gene_labels=True)

markers = ["Pecam1", "Cdh5", "Nos3", "Acta2", "Cnn1", "Tagln", "Rgs5", "Kcnj8", "Col1a1", "Col5a1", "Dpt", "Cd19", "Ighm", "Cd14", "Cd68", "Cd3d"] # Cd3g 없음
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
'MΦ':['Cd14', 'Cd68'],
'T cells':['Cd3d', 'Cd3g']}
reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'M\u03A6', 'T cells')
#reordered = ('Endothelial cells', 'Smooth muscle cells', 'Fibroblasts', 'B cells', 'Macrophages', 'T cells')
test3.obs['celltype'] = test3.obs['celltype'].cat.reorder_categories(list(reordered), ordered=True)

sc.pl.umap(test3, color=['celltype'], add_outline=False, legend_loc='right margin', size=30, color_map='CMRmap', palette='Paired')

# Imputed expression matrix using MAGIC
import magic
test3_MAGIC = test3.copy()
test3_MAGIC.X = test3.layers['scran_log1p']
test3_MAGIC = magic.MAGIC().fit_transform(test3_MAGIC)

test3.layers['magic'] = test3_MAGIC.X
sc.pl.umap(test3, layer='magic', color=['Pecam1'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3, color=['Pecam1'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')


dp = sc.pl.dotplot(test3, layer='magic', celltype_marker, groupby='celltype', return_fig=True)
dp.add_totals(size=1.5, color=['#a6cee3', '#b2df8a', '#fb9a99', '#ff7f00', '#6a3d9a', '#b15928']).legend(colorbar_title='log(SizeFactorNormlized+1)', width=1.5, show_size_legend=False, show_colorbar=False).style(cmap='winter', dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5, grid=True, x_padding=0.4, y_padding=0.6).swap_axes().show()
# MAGIC imputed expression
dp = sc.pl.dotplot(test3, celltype_marker, layer='magic', groupby='celltype', return_fig=True)
dp.add_totals(size=1.5, color=['#a6cee3', '#b2df8a', '#fb9a99', '#ff7f00', '#6a3d9a', '#b15928']).legend(colorbar_title='log(SizeFactorNormlized+1)', width=1.5, show_size_legend=False, show_colorbar=False).style(cmap='winter', dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5, grid=True, x_padding=0.4, y_padding=0.6).swap_axes().show()
sc.pl.StackedViolin(test3, celltype_marker, layer='magic', groupby='celltype')
#test3.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")
#test3 = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3.h5ad")

# Only Endothelial cells
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

test3_endo = anndata.AnnData(X=test3[test3.obs['leiden_r05'].isin(['5', '7'])].layers['counts'], obs=test3[test3.obs['leiden_r05'].isin(['5', '7'])].obs, var=test3[test3.obs['leiden_r05'].isin(['5', '7'])].var)
test3_endo.layers["counts"] = test3_endo.X.copy()

adata_pp = test3_endo.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_endo.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_endo.obs['size_factors'] = size_factors

test3_endo.X /= test3_endo.obs['size_factors'].values[:, None]
test3_endo.X = scipy.sparse.csr_matrix(test3_endo.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)

test3_endo.layers['scran'] = test3_endo.X

sc.pp.log1p(test3_endo) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_endo.layers['scran_log1p'] = test3_endo.X

test3_endo.raw = test3_endo ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

sc.pp.highly_variable_genes(test3_endo)

test3_endo.var['highly_variable'].value_counts() # 2,612 ==> 2021-08-20, # 2,941 ==> 2021-09-28

sc.pp.filter_genes(test3_endo, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE

sc.pp.scale(test3_endo, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
# adata.raw.X의 mean 과 std를 output함
sc.tl.pca(test3_endo, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_endo, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None) #####
sc.tl.umap(test3_endo, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
#test3_endo.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_endo, resolution=0.5, key_added='endo_leiden_r05')
sc.tl.leiden(test3_endo, resolution=1.0, key_added='endo_leiden_r10')

test3_endo.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']

sc.pl.umap(test3_endo, color=['endo_leiden_r05', 'endo_leiden_r10', 'leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3_endo, color=['batch', 'phase', 'percent_mito'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.umap(test3_endo, color=['batch'], group_by='Month1', add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_rank_genes_groups')
#sc.pl.rank_genes_groups(test3_endo, n_genes=5, sharey=False)
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, min_logfoldchange=2, cmap='cividis', show_gene_labels=True, key='endo_leiden_r05_rank_genes_groups')

test3_endo.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

markers = ["Pecam1", "Cdh5", "Nos3", "Acta2", "Cnn1", "Tagln", "Rgs5", "Kcnj8", "Col1a1", "Col5a1", "Dpt", "Cd19", "Ighm", "Cd14", "Cd68", "Cd3d"] # Cd3g 없음
sc.pl.stacked_violin(test4, markers, groupby='batch')

result = test3_endo.uns['endo_leiden_r05_rank_genes_groups']
groups = result['names'].dtype.names
#pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']}).head(5)
#deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pts', 'pts_rest', 'pvals_adj']})
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})
deg_wilcoxon.to_csv("20210930_scanpy_deg_endo.csv", mode='w')

lin = ('0','1','2','3','4','5')
lin = ('1', '3', '5', '4', '2', '0')
#test3_endo.obs['leiden_r05']
test3_endo.obs['endo_leiden_r05'] = test3_endo.obs['endo_leiden_r05'].cat.reorder_categories(list(lin), ordered=True)
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, key='rank_genes_groups', show_gene_labels=True, min_logfoldchange=1, dendrogram=False)

sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_rank_genes_groups')
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, key='endo_leiden_r05_rank_genes_groups', show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap='cividis')
# 아래와 같이 groups을 지정하면 이 group의 순서대로 column이 배치됨 !!!!!!!!!!!!!!!!!!!!!!!! ==> 아닌 것 같다 다시 한 번 확인.
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, key='endo_leiden_r05_rank_genes_groups', groups=['1', '3', '5', '4', '2', '0'], show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap='cividis')



ec_litvinukova = ['Rgcc', 'Car4', 'Sema3g', 'Gja5', 'Plvap']
sc.pl.heatmap(test3_endo, var_names=ec_litvinukova, groupby='leiden_r05', cmap='coolwarm')

artery_schupp = ['Ltbp4', 'Fbln5', 'Bmx', 'Gja4', 'Gja5', 'Efnb2', 'Sox17', 'Sema3g', 'Hey1']
sc.pl.heatmap(test3_endo, var_names=artery_schupp, groupby='leiden_r05', cmap='coolwarm')

capillary_schupp = ['Car4', 'Rgcc', 'Sgk1', 'Sparc', 'Prx']
sc.pl.heatmap(test3_endo, var_names=capillary_schupp, groupby='leiden_r05', cmap='coolwarm')

vein_schupp = ['Nr2f2', 'Selp', 'Vcam1']
sc.pl.heatmap(test3_endo, var_names=vein_schupp, groupby='leiden_r05', cmap='coolwarm')


a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['leiden_r05'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_leiden_r05'] = c
sc.tl.rank_genes_groups(test3_endo, 'aging_leiden_r05', method='wilcoxon', pts=True, key_added='rank_genes_groups_batch_aging')
# ... storing 'aging_leiden_r05' as categorical ==> 만들어 놓은 aging_leiden_r05를 categorical하게 바꿔야함

lin = ('m01_0', 'm10_0', 'm20_0', 'm01_1', 'm10_1', 'm20_1', 'm01_2', 'm10_2', 'm20_2', 'm01_3', 'm10_3', 'm20_3', 'm01_4', 'm10_4', 'm20_4', 'm01_5', 'm10_5', 'm20_5')
#lin = ('m01_0', 'm01_1', 'm01_2', 'm01_3', 'm01_4','m01_5', 'm10_0', 'm10_1', 'm10_2', 'm10_3', 'm10_4','m10_5', 'm20_0', 'm20_1', 'm20_2', 'm20_3', 'm20_4','m20_5')
test3_endo.obs['aging_leiden_r05']
test3_endo.obs['aging_leiden_r05'] = test3_endo.obs['aging_leiden_r05'].cat.reorder_categories(list(lin), ordered=True)

sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=10, groups=['m01_0', 'm10_0', 'm20_0'], groupby='aging_leiden_r05', key='rank_genes_groups_batch_aging', show_gene_labels=True, min_logfoldchange=1, dendrogram=False, cmap='coolwarm')

sc.tl.rank_genes_groups(test3_endo, 'batch', method='wilcoxon', pts=True, key_added='rank_genes_groups_batch')
sc.pl.rank_genes_groups_heatmap(test3_endo, n_genes=15, key='rank_genes_groups_batch', show_gene_labels=True, min_logfoldchange=1)


# Diffusion pseudotime
sc.tl.diffmap(test3_endo)
sc.pl.diffmap(test3_endo, color=['batch', 'Pecam1', 'Cdh5'], add_outline=False, legend_loc='right margin', size=70, color_map='CMRmap')

sc.tl.draw_graph(test3, layout='fa', init_pos=None, neighbors_key=None) ## init_pos가 .obsm에 있는 pca, umap, paga 등이 될 수 있다.
sc.pl.draw_graph(test3, color=['batch', 'PECAM1', 'CDH5', 'phase'], add_outline=True, legend_loc='right margin', size=10, color_map='CMRmap')

start_cell = np.isin(test3_endo.obs['endo_leiden_r05'], '0') # boolean numpy array ==> array([False, False, False, ..., False, False, False])
#max_start_id = np.argmin(test3_endo.obsm['X_diffmap'][start_cell,1]) # 262
max_start_id = np.argmax(test3_endo.obsm['X_diffmap'][start_cell,1])
root_id = np.arange(len(start_cell))[start_cell][max_start_id] # 341
test3_endo.uns['iroot'] = root_id

sc.tl.dpt(test3_endo, n_branchings=1, n_dcs=10) # n_branchings를 0으로 하면 (recommended by Scanpy developer) dpt_groups가 생성 안 됨.
#computing Diffusion Pseudotime using n_dcs=10
sc.pl.dpt_groups_pseudotime(test3_endo) # 여기에서 pseudotime trajecgory 확인.

lin = ('2', '0', '3', '1') # DPT pseudotime group ordering에 맞게 배치
test3_endo.obs['dpt_groups'] = test3_endo.obs['dpt_groups'].cat.reorder_categories(list(lin), ordered=True)
sc.pl.dpt_groups_pseudotime(test3_endo) # 다시 ordering에 맞게 plotting
sc.pl.dpt_timeseries(test3_endo[:, test3_endo.var.highly_variable])


################## aEC 에서 비교 ##################
# endo_leiden_r05에서 '0'과 '2'의 비교? '2'는 m01과 m20이 enrich를 판단하기 어렵고, '0'은 확실히 m01이 많음.

# 0 vs 2 (upregulated in 0)
sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_0vs2_rank_genes_groups', groups=['0'], reference='2')
sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])], n_genes=20, groups=['0'], key='endo_leiden_r05_0vs2_rank_genes_groups', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis')

# 2 vs 0 (upregulated in 2)
sc.tl.rank_genes_groups(test3_endo, 'endo_leiden_r05', method='wilcoxon', pts=True, key_added='endo_leiden_r05_2vs0_rank_genes_groups', groups=['2'], reference='0')
sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])], n_genes=20, groups=['2'], key='endo_leiden_r05_2vs0_rank_genes_groups', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis')












# Only "Arterial" Endothelial cells
test3_aEC = anndata.AnnData(X=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])].layers['counts'], obs=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])].obs, var=test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '2'])].var)
test3_aEC.layers["counts"] = test3_aEC.X.copy()

adata_pp = test3_aEC.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_aEC.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))




del adata_pp
test3_aEC.obs['size_factors'] = size_factors

test3_aEC.X /= test3_aEC.obs['size_factors'].values[:, None]
test3_aEC.X = scipy.sparse.csr_matrix(test3_aEC.X) #왜 이게 새로 들어가야될까?????
test3_aEC.X

sc.pp.log1p(test3_aEC) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_aEC.raw = test3_aEC ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

###################################### 2. BBKNN  ######################################
sc.pp.highly_variable_genes(test3_aEC)

test3_aEC.var['highly_variable'].value_counts() # 2,798 ==> 2021-10-14

sc.pp.scale(test3_aEC, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
sc.tl.pca(test3_aEC, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_aEC, batch_key='batch', n_pcs=7, neighbors_within_batch=3, trim=None)
sc.tl.umap(test3_aEC, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_aEC.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_aEC, resolution=0.5, key_added='aEC_leiden_r05')
sc.pl.umap(test3_aEC, color=['batch', 'aEC_leiden_r05', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

#test3_aEC.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")
#test3_aEC = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")


# Diffusion pseudotime (########################### Fail ###########################)
sc.tl.diffmap(test3_aEC, n_comps=50) # n_comps=15 (default)
sc.pl.diffmap(test3_aEC, color=['batch', 'Pecam1', 'Cdh5'], add_outline=False, legend_loc='right margin', size=70, color_map='CMRmap')

start_cell = np.isin(test3_aEC.obs['aEC_leiden_r05'], '0') # boolean numpy array ==> array([False, False, False, ..., False, False, False])
max_start_id = np.argmax(test3_aEC.obsm['X_diffmap'][start_cell,2])
root_id = np.arange(len(start_cell))[start_cell][max_start_id] # 253
test3_aEC.uns['iroot'] = root_id

sc.tl.dpt(test3_aEC, n_branchings=1, n_dcs=50) # n_branchings를 0으로 하면 (recommended by Scanpy developer) dpt_groups가 생성 안 됨.
#computing Diffusion Pseudotime using n_dcs=10
sc.pl.dpt_groups_pseudotime(test3_aEC) # 여기에서 pseudotime trajecgory 확인.

lin = ('1', '0', '2') # DPT pseudotime group ordering에 맞게 배치
test3_aEC.obs['dpt_groups'] = test3_aEC.obs['dpt_groups'].cat.reorder_categories(list(lin), ordered=True)
sc.pl.dpt_groups_pseudotime(test3_aEC) # 다시 ordering에 맞게 plotting
sc.pl.dpt_timeseries(test3_aEC[:, test3_aEC.var.highly_variable])

# PAGA

sc.tl.paga(test3_aEC, groups='aEC_leiden_r05')
sc.pl.paga(test3_aEC, color=['aEC_leiden_r05'], threshold=0.2)
#--> added 'pos', the PAGA positions (adata.uns['paga'])
sc.tl.umap(test3_aEC, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='paga', method='umap')
sc.pl.umap(test3_aEC, color=['batch', 'aEC_leiden_r05', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')
sc.pl.paga_compare(test3_aEC, threshold=0.02, size=10, frameon=True, edges=True)











# age difference ####################################################
a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['endo_leiden_r05'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_endo_leiden_r05'] = c
sc.tl.rank_genes_groups(test3_endo, 'aging_endo_leiden_r05', method='wilcoxon', groups=['m01_0', 'm10_0', 'm20_0'], pts=True, key_added='rank_genes_groups_aging0')
sc.tl.rank_genes_groups(test3_endo, 'aging_endo_leiden_r05', method='wilcoxon', pts=True, key_added='rank_genes_groups_aging2')


sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['aging_endo_leiden_r05'].isin(['m01_0', 'm10_0', 'm20_0'])], n_genes=20, groups=['m01_0', 'm10_0', 'm20_0'], key='rank_genes_groups_aging0', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis', use_raw=False)







# Removing vasa vasorum, lymphatic EC and SMC-like EC
test3_aEC = anndata.AnnData(X=test3_endo[~test3_endo.obs['endo_leiden_r05'].isin(['4', '5'])].layers['counts'], obs=test3_endo[~test3_endo.obs['endo_leiden_r05'].isin(['4', '5'])].obs, var=test3_endo[~test3_endo.obs['endo_leiden_r05'].isin(['4', '5'])].var)
test3_aEC.layers["counts"] = test3_aEC.X.copy()

adata_pp = test3_aEC.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_aEC.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))




del adata_pp
test3_aEC.obs['size_factors'] = size_factors

test3_aEC.X /= test3_aEC.obs['size_factors'].values[:, None]
test3_aEC.X = scipy.sparse.csr_matrix(test3_aEC.X) #왜 이게 새로 들어가야될까?????
test3_aEC.X

sc.pp.log1p(test3_aEC) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_aEC.raw = test3_aEC ## ==> log transforamtion 된 것이 raw로 들어가게 됨.

###################################### 2. BBKNN  ######################################
sc.pp.highly_variable_genes(test3_aEC)

test3_aEC.var['highly_variable'].value_counts() # 2,798 ==> 2021-10-14

sc.pp.scale(test3_aEC, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
sc.tl.pca(test3_aEC, n_comps=100, use_highly_variable=True, svd_solver='arpack')

#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_aEC, batch_key='batch', n_pcs=7, neighbors_within_batch=3, trim=None)
sc.tl.umap(test3_aEC, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_aEC.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_aEC, resolution=0.5, key_added='aEC_leiden_r05')
sc.pl.umap(test3_aEC, color=['batch', 'aEC_leiden_r05', 'endo_leiden_r05'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')

#test3_aEC.write(filename="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")
#test3_aEC = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_aEC.h5ad")







############# cellrank (CytoTrace) ######################

# on test3_endo
import scvelo as scv
from cellrank.tl.kernels import CytoTRACEKernel

# scVelo hack
test3_endo.layers['spliced'] = test3_endo.raw.X
test3_endo.layers['unspliced'] = test3_endo.raw.X

scv.pp.moments(test3_endo)
ctk_endo = CytoTRACEKernel(test3_endo)
sc.pl.umap(test3_endo, color=['ct_pseudotime'], add_outline=False, legend_loc='right margin', size=150)

ctk_endo.compute_transition_matrix(threshold_scheme="soft", nu=0.5)
ctk_endo.compute_projection(basis="umap")

scv.pl.velocity_embedding_stream(test3_endo, color="batch", vkey="T_fwd", basis="umap", legend_loc="right")

from cellrank.tl.estimators import GPCCA
g_fwd = GPCCA(ctk_endo)
g_fwd.compute_schur(n_components=20)
g_fwd.plot_spectrum(real_only=True)









############################################## cellrank ##############################################################

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
plt.rcParams['figure.figsize'] = (7,7)

import cellrank as cr
from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA

test3_endo = sc.read("test3_endo.h5ad")
timepoint = [1] * batches.count('m01') + [10] * batches.count('m10') + [20] * batches.count('m20')
test3_endo.obs['months'] = timepoint

wk = WOTKernel(test3_endo, time_key="months")

wk.compute_initial_growth_rates(organism="mouse", key_added="growth_rate_init")
#WARNING: genes are not in var_names and ignored: ['Hn1', 'Mlf1ip', 'Fam64a']
#WARNING: genes are not in var_names and ignored: ['Adck3', 'Nhlh2', 'Ikbkap', 'Gnb2l1', 'Krt17']

wk.compute_transition_matrix(growth_iters=3, growth_rate_key="growth_rate_init", last_time_point="connectivities")
# ==> ERROR: TypeError: compute_transport_map() got an unexpected keyword argument 'cost_matrix'









############################################## MAGIC ##############################################################

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

test3_endo = anndata.AnnData(X=test3[test3.obs['leiden_r05'].isin(['5', '7'])].layers['counts'], obs=test3[test3.obs['leiden_r05'].isin(['5', '7'])].obs, var=test3[test3.obs['leiden_r05'].isin(['5', '7'])].var)
test3_endo.layers["counts"] = test3_endo.X.copy()

adata_pp = test3_endo.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.tl.pca(adata_pp, n_comps=15) ## 여기서 이 n_component의 숫자를 늘리면 size_factors를 estimation하는 데 도움이 될까?
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)
input_groups = adata_pp.obs['groups']
data_mat = test3_endo.X.T
%%R -i data_mat -i input_groups -o size_factors
size_factors = BiocGenerics::sizeFactors(computeSumFactors(SingleCellExperiment::SingleCellExperiment(list(counts=data_mat)), clusters=input_groups, min.mean=0.1))



del adata_pp
test3_endo.obs['size_factors'] = size_factors
test3_endo.X /= test3_endo.obs['size_factors'].values[:, None]
test3_endo.X = scipy.sparse.csr_matrix(test3_endo.X) #왜 이게 새로 들어가야될까????? # 아니면 ERRROR 남 (highly_variable_genes에서)
test3_endo.X
sc.pp.log1p(test3_endo) # works on anndata.X
#integrated.X = scipy.sparse.csr_matrix(integrated.X)
test3_endo.raw = test3_endo ## ==> log transforamtion 된 것이 raw로 들어가게 됨.
sc.pp.highly_variable_genes(test3_endo)
test3_endo.var['highly_variable'].value_counts() # 2,612 ==> 2021-08-20, # 2,941 ==> 2021-09-28
sc.pp.filter_genes(test3_endo, min_cells=0) # integrated.var에 n_cells 추가 ==> test3에서 이루어졌던 n_cells UPDATE
sc.pp.scale(test3_endo, max_value=10) # ... as `zero_center=True`, sparse input is densified and may lead to large memory consumption
sc.tl.pca(test3_endo, n_comps=100, use_highly_variable=True, svd_solver='arpack')
#sce.pp.bbknn default ==> n_pcs=50, neighbors_within_batch=3, trim=None, annoy_n_trees=10,
sce.pp.bbknn(test3_endo, batch_key='batch', n_pcs=20, neighbors_within_batch=5, trim=None) #####
sc.tl.umap(test3_endo, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
test3_endo.uns['batch_colors'] = ['#2a2b2d', '#2da8d8', '#d9514e']
sc.tl.leiden(test3_endo, resolution=0.5, key_added='endo_leiden_r05')
test3_endo.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']

import magic
test3_endo_MAGIC = test3_endo.copy()
test3_endo_MAGIC.X = test3_endo_MAGIC.layers['scran_log1p']
test3_endo_MAGIC = magic.MAGIC().fit_transform(test3_endo_MAGIC)
test3_endo.layers['magic'] = test3_endo_MAGIC.X
del test3_endo_MAGIC


test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")
test3_endo.layers['magic'] = test3_endo_MAGIC.X
sc.pl.umap(test3_endo, layer='magic', color=['Pecam1'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')



############################################## Slingshot ##############################################################
library(slingshot)
library(tidyverse)
library(tidymodels)

sds <- slingshot(Embeddings(seu.subset.endo, "umap"), clusterLabels=seu.subset.endo$seurat_clusters, start.clus=3, stretch=2)
# A Function (that assigns colors to each cell in base R graphics https://bustools.github.io/BUS_notebooks_R/slingshot.html)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

library(scales)
cell_colors <- cell_pal(seu.subset.endo$orig.ident, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu.subset.endo$seurat_clusters, hue_pal())
plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 1)
#lines(sds, lwd = 2, type = 'lineages', col = 'black')
lines(sds, lwd = 2, col = 'black') ###

library(viridis)
nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end=0.95)

##
par(mfrow = c(nr, nc))
for (i in nms) {
	colors <- pal[cut(pt[,i], breaks=100)]
	plot(reducedDim(sds), col=colors, pch=16, cex=1, main=i)
	#lines(sds, lwd=2, col='black', type='lineages')
	lines(sds, lwd=2, col='black')
}
##





top_hvg <- VariableFeatures(seu.subset.endo)
dat_use <- t(GetAssayData(seu.subset.endo, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use)
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>% set_engine("ranger", importance = "impurity", num.threads = 3) %>% fit(pseudotime ~ ., data = dat_train)
val_results <- dat_val %>% mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% select(truth = pseudotime, estimate)
metrics(data = val_results, truth, estimate)








#### Table generation 2021-09-16

a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['endo_leiden_r05'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_endo_leiden_r05'] = c
lin = ('m01_0', 'm10_0', 'm20_0', 'm01_1', 'm10_1', 'm20_1', 'm01_2', 'm10_2', 'm20_2', 'm01_3', 'm10_3', 'm20_3', 'm01_4', 'm10_4', 'm20_4', 'm01_5', 'm10_5', 'm20_5')
#lin = ('m01_0', 'm01_1', 'm01_2', 'm01_3', 'm01_4','m01_5', 'm10_0', 'm10_1', 'm10_2', 'm10_3', 'm10_4','m10_5', 'm20_0', 'm20_1', 'm20_2', 'm20_3', 'm20_4','m20_5')
test3_endo.obs['aging_endo_leiden_r05']
test3_endo.obs['aging_endo_leiden_r05'] = test3_endo.obs['aging_endo_leiden_r05'].astype('category').cat.reorder_categories(list(lin), ordered=True)

tipcell_markers = ['Kdr', 'Flt4', 'Nrp1', 'Nrp2', 'Pdgfb', 'Dll4', 'Angpt2', 'Apln', 'Unc5b', 'Robo4', 'Plxnd1', 'Efnb2', 'Cxcr4']
sc.tl.score_genes(test3_endo, tipcell_markers, score_name='tipcell_score', use_raw=True)
sc.pl.umap(test3_endo, layer='magic', color=tipcell_markers, color_map=cmap, ncols=3)

