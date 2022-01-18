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


integrated = AnnData.concatenate(m03_1, m03_2, m03_3, m03_4, join='outer', batch_categories = ['m01', 'm10', 'm20'], index_unique = '-')
