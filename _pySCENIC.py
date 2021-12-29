#### pySCENIC
from anndata import AnnData
import anndata
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sb
import scanpy as sc
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (7,7)

matplotlib.use('TkAgg')
%matplotlib

# STEP 1: Gene regulatory network inference, and generation of co-expression modules
import loompy as lp
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")
row_attrs = {"Gene": np.array(test3_endo.var_names) ,}
col_attrs = {"CellID": np.array(test3_endo.obs_names),
             "nGene": np.array(np.sum(test3_endo.X.transpose()>0 , axis=0)).flatten(),
             "nUMI": np.array(np.sum(test3_endo.X.transpose(), axis=0)).flatten(),}
lp.create("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/test3_endo.loom", test3_endo.X.transpose(), row_attrs, col_attrs)

f_loom_path_scenic="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/test3_endo.loom"
f_tfs="/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/pySCENIC_Database/mm_mgi_tfs.txt"

!pyscenic grn \
{f_loom_path_scenic} \
{f_tfs} \
-o test3_endo_adj.csv \
--num_workers 20

# STEP 2-3: Regulon prediction aka cisTarget from CLI
!pyscenic ctx test3_endo_adj.csv \
/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname /mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
--expression_mtx_fname {f_loom_path_scenic} \
--output reg.csv \
--mask_dropouts \
--num_workers 20

# STEP 4: Cellular enrichment (aka AUCell) from CLI

nGenesDetectedPerCell = np.sum(test3_endo.X>0, axis=1)
# test3_endo.X는 np.shape(test3_endo.X) ==> (879,19944) 즉, row=cell, column=gene으로 구성되어 있는데,
# 여기에서 양의 값을 가지는 element만 골라서 axis=1 에서 더한다.
# 주의할것은 axis=0은 2D array에서 각 열에 해당하는 행을 의미하고, axis=1은 각 행에 해당하는 열을 의미함.
# 즉, a = np.array([[1,2,3],[4,5,6]]) 에서 np.sum(a, axis=0) ==> array([5, 7, 9])이고, np.sum(a, axis=1) ==> array([ 6, 15]) 이다.
# 즉 위의 nGenesDetectedPerCell에서는 element가 양수값을 가지는 것들만 모아서 각각의 cell에서 detected 된 유전자 갯수를 셀 수 있다.
sb.histplot(nGenesDetectedPerCell, bins='fd')

percentiles = np.quantile(nGenesDetectedPerCell, [.01, .05, .10, .50, 1])
# 여기에서 유전자 갯수들과 그것의 histogram이 있을텐데, np.quantile을 이용하여 histogram의 y value 1퍼센트 5퍼센트 등에 해당하는 유전자 값을 구할 수 있을 것이다.
percentiles
#array([1440.9, 1558.4, 1679. , 2815. , 6307. ])
#즉 1558.4개의 유전자가 cell의 5%를 차지한다는 것을 알 수 있다. ==> 이게 alxdml aucell의 5% default value가 됨

!pyscenic aucell \
    test3_endo.loom \
    test3_endo_reg.csv \
    --num_workers 30 \
    --output test3_endo_pyscenic_output.loom

import loompy as lp
test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")

lf = lp.connect("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC/test3_endo_pyscenic_output.loom", mode='r+', validate=False)
lf.ca.keys()
#['CellID', 'RegulonsAUC', 'nGene', 'nUMI']
lf.ra.keys()
#['Gene', 'Regulons']
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)


import umap
import json
import zlib
import base64

runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform # parameter 조정! ==> 이거 꼭 필요할 것 같음!!
dr_umap = runUmap(auc_mtx)
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv("scenic_umap.txt", sep='\t')

# Integrate output the result from pySCENIC and the scanpy analysis
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
regulons = lf.ra.Regulons
dr_umap = pd.read_csv('/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC/scenic_umap.txt', sep='\t', header=0, index_col=0)

auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(') # columns: 2010315B03Rik(+) ==> 2010315B03Rik_(+)
regulons.dtype.names = tuple([ x.replace("(","_(") for x in regulons.dtype.names])
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
tmp = x.get('regulon').replace("(","_(")
x.update({'regulon': tmp})


umapDF = pd.DataFrame(test3_endo.obsm['X_umap'], columns=['_X', '_Y'])
Embeddings_X = pd.DataFrame(index=lf.ca.CellID)
Embeddings_X = pd.concat([pd.DataFrame(test3_endo.obsm['X_umap'], index=test3_endo.obs.index)[0],
                          pd.DataFrame(test3_endo.obsm['X_pca'], index=test3_endo.obs.index)[0],
                          dr_umap['X']],
                          sort=False, axis=1, join='outer')
Embeddings_X.columns = ['1','2','3'] # '1': scanpy UMAP coordinate, '2': scanpy PCA coordinate, '3': pyscenic UMAP

Embeddings_Y = pd.DataFrame(index=lf.ca.CellID)
Embeddings_Y = pd.concat([pd.DataFrame(test3_endo.obsm['X_umap'], index=test3_endo.obs.index)[1],
                          pd.DataFrame(test3_endo.obsm['X_pca'], index=test3_endo.obs.index)[1],
                          dr_umap['Y']],
                          sort=False, axis=1, join='outer')
Embeddings_Y.columns = ['1','2','3'] # '1': scanpy UMAP coordinate, '2': scanpy PCA coordinate, '3': pyscenic UMAP

metaJson = dict()
metaJson['embeddings'] = [{"id": 1, "name": f"Scanpy UMAP (HVGs)"},
                          {"id": 2, "name": "Scanpy PC1/PC2"},
                          {"id": 3, "name": "pySCENIC AUC UMAP"}]
metaJson['clusterings'] = [{"id": 0, "group": "Scanpy", "name": "Scanpy leiden resolution 0.5", "clusters":[]}]
metaJson['metrics'] = [{"name": "nUMI"},
                       {"name": "nGene"},]
# {"name": "percent_mito"}
metaJson['annotations'] = [{"name": "Leiden_clusters", "values": list(set(test3_endo.obs['endo_leiden_r05'].astype(str)))},
                           {"name": "Timepoint", "values": list(set(test3_endo.obs['batch'].values))},]

# pySCENIC regulon thresholds
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in test3_endo.obs['endo_leiden_r05']])) + 1): # range(0,6) ==> 0,1,2,3,4,5 (test3_endo end_leiden_r05)
clustDict = dict()
clustDict['id'] = i
clustDict['description'] = f'Unannotated Cluster {i+1}' # 얘 i+1이 아니라 i라고 해야할 듯
metaJson['clusterings'][0]['clusters'].append(clustDict)



clusterings = pd.DataFrame()
clusterings["0"] = test3_endo.obs['endo_leiden_r05'].values.astype(np.int64)

def dfToNamedMatrix(df):
array_pre = [tuple(i) for i in df.values]
dtype = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
array = np.array(array_pre, dtype=dtype)
return array

col_attrs = dict(CellID = np.array(test3_endo.obs.index),
                 nUMI = np.array(test3_endo.obs['n_counts'].values),
                 nGene = np.array(test3_endo.obs['n_genes'].values),
                 Leiden_clusters = np.array(test3_endo.obs['endo_leiden_r05'].values),
                 Timepoint = np.array(test3_endo.obs['batch'].values),
                 Embedding = dfToNamedMatrix(umapDF),
                 Embeddings_X = dfToNamedMatrix(Embeddings_X),
                 Embeddings_Y = dfToNamedMatrix(Embeddings_Y),
                 RegulonsAUC = dfToNamedMatrix(auc_mtx),
                 Clusterings = dfToNamedMatrix(clusterings),
                 ClusterID = np.array(test3_endo.obs['endo_leiden_r05'].values))

row_attrs = dict(Gene = lf.ra.Gene,
                 Regulons = regulons)

attrs = dict(title = "sampleTitle",
             MetaData = json.dumps(metaJson),
             Genome = 'mm10',
             SCopeTreeL1 = "",
             SCopeTreeL2 = "",
             SCopeTreeL3 = "")
# json.dumps(metaJson) ==> metaJson (type: dict)을 str으로

# compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

# Create a new loom file, copying the expression matrix from the open loom connection:
lp.create(filename="test3_endo_pyscenic_integrated_output.loom",
          layers=lf[:,:],
          row_attrs=row_attrs,
          col_attrs=col_attrs,
          file_attrs=attrs)
lf.close()

# Clustering the data
sb.clustermap(auc_mtx, method='ward', metric='euclidean', z_score=None, cmap="rocket_r")
sb.clustermap(auc_mtx, method='ward', metric='euclidean', z_score=None, standard_scale=1, cmap="rocket_r", cbar_pos=(0,0,0,0))


#### Downstream analysis ####

#from pyscenic.plotting import plot_binarization
#from pyscenic.export import add_scenic_metadata
#from pyscenic.cli.utils import load_signatures

# pySCENIC output
lf = lp.connect("test3_endo_pyscenic_integrated_output.loom", mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData)))
exprMat = pd.DataFrame(lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

# Creating a dictionary of regulons
regulons = dict()
for i,r in pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene).iteritems():
# i ==> 이 dataframe의 column name (Regulon name), r ==> gene name을 row index로 하는 각각의 regulon value의 pandas series
regulons[i] = list(r[r==1].index.values)

# cell annotations from the loom column attributes
lf.ca.keys()
# ['CellID', 'ClusterID','Clusterings','Embedding','Embeddings_X','Embeddings_Y','Leiden_clusters','RegulonsAUC','Timepoint','nGene', 'nUMI']
cellAnnot = pd.concat([pd.DataFrame(lf.ca.ClusterID, index=lf.ca.CellID),
                      pd.DataFrame(lf.ca.Leiden_clusters, index=lf.ca.CellID),
                      pd.DataFrame(lf.ca.Timepoint, index=lf.ca.CellID),
                      pd.DataFrame(lf.ca.nGene, index=lf.ca.CellID),
                      pd.DataFrame(lf.ca.nUMI, index=lf.ca.CellID)], axis=1)
cellAnnot.columns = ['ClusterID', 'Leiden_clusters', 'Timepoint', 'nGene', 'nUMI']

# Capturing Embeddings
dr = [pd.DataFrame(lf.ca.Embedding, index=lf.ca.CellID)] # 얘가 redundant함. drx와 dry에서 UMAP 두번 처리됨 ==> 확인해보셈 진짜임
dr_names = [meta['embeddings'][0]['name'].replace(" ", "_")]
# meta['embeddings']
# [{'id': 1, 'name': 'Scanpy UMAP (HVGs)'},
#  {'id': 2, 'name': 'Scanpy PC1/PC2'},
#  {'id': 3, 'name': 'pySCENIC AUC UMAP'}]
drx = pd.DataFrame(lf.ca.Embeddings_X, index=lf.ca.CellID)
dry = pd.DataFrame(lf.ca.Embeddings_Y, index=lf.ca.CellID)

for i in range( len(drx.columns) ):
dr.append(pd.concat([drx.iloc[:, i], dry.iloc[:, i]], sort=False, axis=1, join='outer'))
dr_names.append( meta['embeddings'][i]['name'].replace(" ","_").replace('/','-') )

# Rename the dimensional reduction (dr) columns
for i, x in enumerate(dr):
x.columns = ['X', 'Y']

for i,x in enumerate(dr):
    test3_endo.obsm['X_' + dr_names[i]] = x.to_numpy()

#test3_endo.obsm
#AxisArrays with keys: X_pca, X_umap, X_Scanpy_UMAP_(HVGs), X_Scanpy_PC1-PC2, X_pySCENIC_AUC_UMAP

test3_endo.uns['batch_colors'] = ['#689aff', '#fdbf6f', '#b15928']
sc.pl.scatter(test3_endo, basis='pySCENIC_AUC_UMAP', color=['endo_leiden_r05', 'batch'])


### Regulon specificity scores (RSS) across endo leiden clusters (resolution = 0.5)

cellAnnot.to_csv("test3_endo_cellAnnot.txt", sep='\t')
auc_mtx.to_csv("test3_endo_auc_mtx.txt", sep='\t')
# source activate pySCENIC

from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
import pandas as pd

cellAnnot = pd.read_csv('test3_endo_cellAnnot.txt', sep='\t', header=0, index_col=0)
auc_mtx = pd.read_csv('test3_endo_auc_mtx.txt', sep='\t', header=0, index_col=0)
# cellAnnot column names ==> ClusterID  Leiden_clusters Timepoint  nGene     nUMI

# Calculating Regulon specificity scores (RSS) acorss predicted cell types
rss_cellType = regulon_specificity_scores(auc_mtx, cellAnnot['Leiden_clusters'])

# RSS panel plot
cats = sorted(list(set(cellAnnot['Leiden_clusters'])))
fig = plt.figure(figsize=(15, 8))
for c, num in zip(cats, range(1, len(cats)+1)):
x = rss_cellType.T[c]
ax = fig.add_subplot(1,5,num) # 1X5 grid의 num 번째에 subplot을 추가
plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05, x.max() + (x.max() - x.min()) * 0.05)
for t in ax.texts:
t.set_fontsize(10)
adjust_text(ax.texts, autoalign='xy', ha='right', va='top', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001)


# ha ==> {'center', 'right', 'left'}
# va ==> {'center', 'top', 'bottom', 'baseline', 'center_baseline'}
rss_cellType.T[0].sort_values(ascending=False)

# Select the top 5 regulons from each cell type
topreg = list()
for i,c in enumerate(cats):
topreg.extend(list(rss_cellType.T[c].sort_values(ascending=False)[:5].index))
# Generate a Z-score for each regulon to enable comparison between regulons
auc_mtx_Z = pd.DataFrame(index=auc_mtx.index)
for col in list(auc_mtx.columns):
auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)


# Heatmap generation for topreg
colors1 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'] # leiden cluster
colorsd1 = dict(zip(cats,colors1))
colormap1 = [colorsd1[x] for x in cellAnnot['Leiden_clusters']]

colors2 = ['#689aff', '#fdbf6f', '#b15928'] # Timepoint
colorsd2 = dict(zip(['m01', 'm10', 'm20'], colors2))
colormap2 = [colorsd2[x] for x in cellAnnot['Timepoint']]

colormap3 = list(itertools.chain(*list(map(lambda x: [x]*5, colors1)))) # row colors

sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray', yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=[colormap1, colormap2], col_colors=colormap3, cmap="YlGnBu", figsize=(7,5), method='ward')
g.cax.set_visible(False) # To remove colorbar





# Calculating Regulon specificity scores (RSS) acorss Timepoint
rss_timepoint = regulon_specificity_scores(auc_mtx, cellAnnot['Timepoint'])
# Select the top 5 regulons from each "Timepoint"
topreg2 = list()
cats2 = ['m01', 'm10', 'm20']
for i,c in enumerate(cats2):
topreg2.extend(list(rss_timepoint.T[c].sort_values(ascending=False)[:5].index))
# Generate a Z-score for each regulon to enable comparison between regulons
auc_mtx_Z = pd.DataFrame(index=auc_mtx.index)
for col in list(auc_mtx.columns):
auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

# Heatmap generation for topreg
colors1 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'] # leiden cluster
colorsd1 = dict(zip(cats,colors1))
colormap1 = [colorsd1[x] for x in cellAnnot['Leiden_clusters']]

colors2 = ['#689aff', '#fdbf6f', '#b15928'] # Timepoint
colorsd2 = dict(zip(['m01', 'm10', 'm20'], colors2))
colormap2 = [colorsd2[x] for x in cellAnnot['Timepoint']]

colormap4 = list(itertools.chain(*list(map(lambda x: [x]*5, colors2)))) # row colors

sns.set(font_scale=1.2)
g2 = sns.clustermap(auc_mtx_Z[topreg2], annot=False,  square=False,  linecolor='gray', yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=[colormap1, colormap2], col_colors=colormap4, cmap="YlGnBu", figsize=(7,5))
g2.cax.set_visible(False) # To remove colorbar







from pyscenic.utils import load_motifs
import operator as op
from IPython.display import HTML, display

# View the motifs table
df_motifs = load_motifs("test3_endo_reg.csv")



# Network inference output
import loompy as lp
from pyscenic.utils import modules_from_adjacencies
import pandas as pd

adjacencies = pd.read_csv("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/test3_endo_adj.csv", index_col=False, sep=',')

lf = lp.connect("test3_endo_pyscenic_integrated_output.loom", mode='r', validate=False)
exprMat = pd.DataFrame(lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
regulons = dict()
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

# creating modules
modules = list(modules_from_adjacencies(adjacencies, exprMat))
tf = 'Nfia' # top regulons in test3_endo leiden cluster 0
tf_mods = [x for x in modules if x.transcription_factor==tf] # modules for Transcription Factor (ex. Nfia)

for i, mod in enumerate(tf_mods):
    print(f'{tf} module {str(i)}: {len(mod.genes)} genes')
print(f'{tf} regulon: {len(regulons[tf+"_(+)"])} genes')

for i,mod in enumerate(tf_mods):
    with open(tf+'_module_'+str(i)+'.txt', 'w') as f:
        for item in mod.genes:
            f.write("%s\n" % item)

with open(tf+'_regulon.txt', 'w') as f:
    for item in regulons[tf+'__(+)']:
        f.write("%s\n" % item)
