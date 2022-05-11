##### pySCENIC on EC subclustsers #####


## Environment (scanpy_1.8.1)
from anndata import AnnData
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import loompy as lp
%autoindent

## STEP 1: Gene regulatory network inference, and generation of co-expression modules.
import loompy as lp

test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")
test3_endo = test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])]
# Creating input for STEP 1.
row_attrs = {"Gene": np.array(test3_endo.var_names), }
col_attrs = {"CellID": np.array(test3_endo.obs_names),
             "nGene": np.array(np.sum(test3_endo.layers['counts'].toarray().transpose() > 0, axis=0)).flatten(),
             "nUMI": np.array(np.sum(test3_endo.layers['counts'].toarray().transpose(), axis=0)).flatten(), }

# np.sum(test3_endo.layers['counts'].toarray() > 0, axis=0) ==> 특정 cell 특정 gene에 UMI count가 하나라도 있을 경우 이 gene을 counting (즉, gene 숫자를 세는 것이 됨 equivalent to test3.obs['n_genes'])
# np.sum(test3_endo.layers['counts'].toarray().transpose(), axis=0) ==> 각 cell의 UMI count를 세게 됨 (equivalent to test3_endo.obs['n_counts'])

# 얘는 사실 위에 처럼 저렇게 복잡할 필요가 없음 : .transpose() 안 쓰고 axis=1로 하면 됨.

# row_attrs = {"Gene": np.array(test3_endo.var_names) ,}
# col_attrs = {"CellID": np.array(test3_endo.obs_names),
#             "nGene": np.array(np.sum(test3_endo.layers['counts'].toarray() > 0, axis=1)).flatten(),
#             "nUMI": np.array(np.sum(test3_endo.layers['counts'].toarray() axis=1)).flatten(),}


lp.create("/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo.loom",
          test3_endo.layers['counts'].toarray().transpose(),
          row_attrs,
          col_attrs)

# Environment (pySCENIC)
# pySCENIC Execution # num_workers 20으로 해야만 됨
nohup \
pyscenic grn \
--output test3_endo_adj.csv \
--method grnboost2 \
--seed 42 \
--num_workers 20 \
/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo.loom \
/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/pySCENIC_Database/mm_mgi_tfs.txt > grn_20220221.out &

# Environment (pySCENIC)
## STEP 2 to 3: Regulon prediction aka cisTarget from CLI
nohup \
pyscenic ctx \
test3_endo_adj.csv \
/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname /mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
--expression_mtx_fname /data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo.loom \
--mask_dropouts \
--rank_threshold 5000 \
--auc_threshold 0.05 \
--nes_threshold 3.0 \
--num_workers 20 \
--output test3_endo_reg.csv \
> ctx_20220221.out &

# Environment (scanpy_1.8.1)
## STEP 4: Cellular enrichment (aka AUCell) from CLI
# IMPORTANT to CHECK that most cells have a substantial fraction of expressed/detected genes in the calculation of AUC.
import seaborn as sns
%matplotlib
nGenesDetectedPerCell = np.sum(test3_endo.layers['counts'].toarray() > 0, axis=1)
percentiles = pd.Series(np.percentile(nGenesDetectedPerCell, [1, 5, 10, 50, 100]),
                        index=['1%', '5%', '10%', '50%', '100%'])

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(10, 6))
for i, x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0, ymax=1, color='red')
    ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]})', color='red', rotation=30,
            size='x-small', rotation_mode='anchor')

ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')

# Environment (pySCENIC)
nohup \
pyscenic aucell \
/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo.loom \
/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo_reg.csv \
--rank_threshold 5000 \
--auc_threshold 0.05 \
--nes_threshold 3.0 \
--num_workers 10 \
--seed 42 \
--output test3_endo_pyscenic_output.loom > aucell_20220221.out &

## pySCENIC's AUC matrix retrieval
# Environment (scanpy_1.8.1)
import scanpy as sc
import loompy as lp
import pandas as pd
import numpy as np
import json
import zlib
import base64
%autoindent
%matplotlib
sns.set(font="Arial", font_scale=1, style='ticks')

lf = lp.connect(
    "/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo_pyscenic_output.loom",
    mode='r+', validate=False)
lf.ca.keys()
# ['CellID', 'RegulonsAUC', 'nGene', 'nUMI']
lf.ra.keys()
# ['Gene', 'Regulons']
lf.attrs.keys()
# ['CreationDate', 'LOOM_SPEC_VERSION', 'MetaData', 'last_modified']

auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
# 이 MetaData는 contains regulons and each regulon's threshold values ('defaultThresholdValue'), name for determining threshold ('defaultThresholdName': 'gaussian_mixture_split')
# and lastly, 'motifData' for each regulon as .png file.

# Regulon의 이름들을 다음과 같이 바꿔주는 작업 010315B03Rik(+) ==> 2010315B03Rik_(+) ####
auc_mtx.columns = auc_mtx.columns.str.replace('\(', '_(')  # columns: 2010315B03Rik(+) ==> 2010315B03Rik_(+)
regulons.dtype.names = tuple(
    [x.replace("(", "_(") for x in regulons.dtype.names])  # type( regulons.dtype.names ) == tuple
rt = meta['regulonThresholds']
for i, x in enumerate(rt):
    tmp = x.get('regulon').replace("(", "_(")
    x.update({'regulon': tmp})

# test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")
# test3_endo = test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])]
metaJson = dict()
metaJson['clusterings'] = [{"id": 0, "group": "Scanpy", "name": "Scanpy leiden resolution 0.5", "clusters": []}]
metaJson['metrics'] = [{"name": "nUMI"},
                       {"name": "nGene"}, ]

metaJson['annotations'] = [
    {"name": "Leiden_clusters", "values": list(set(test3_endo.obs['endo_leiden_r05'].astype(str)))},
    {"name": "Timepoint", "values": list(set(test3_endo.obs['batch'].values))}, ]
metaJson["regulonThresholds"] = rt  # pySCENIC regulonThreholds

# 이 아래 for문 문제 있을 수도 있음 나중에 꼭 확인 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for i in range(max(set([int(x) for x in test3_endo.obs[
    'endo_leiden_r05']])) + 1):  # range(0,6) ==> 0,1,2,3,4,5 (test3_endo end_leiden_r05)
    clustDict = dict()
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i}'  # 얘 i+1이 아니라 i라고 해야할 듯
    metaJson['clusterings'][0]['clusters'].append(clustDict)


def dfToNamedMatrix(df):
    array_pre = [tuple(i) for i in df.values]
    dtype = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    array = np.array(array_pre, dtype=dtype)
    return array

clusterings = pd.DataFrame()
clusterings["0"] = test3_endo.obs['endo_leiden_r05'].values.astype(np.int64)

col_attrs = dict(CellID=np.array(test3_endo.obs.index),
                 nUMI=np.array(test3_endo.obs['n_counts'].values),
                 nGene=np.array(test3_endo.obs['n_genes'].values),
                 Leiden_clusters=np.array(test3_endo.obs['endo_leiden_r05'].values),
                 Timepoint=np.array(test3_endo.obs['batch'].values),
                 RegulonsAUC=dfToNamedMatrix(auc_mtx),
                 Clusterings=dfToNamedMatrix(clusterings),
                 ClusterID=np.array(test3_endo.obs['endo_leiden_r05'].values))

row_attrs = dict(Gene=lf.ra.Gene,
                 Regulons=regulons)

attrs = dict(title="sampleTitle",
             MetaData=json.dumps(metaJson),
             Genome='mm10',
             SCopeTreeL1="",
             SCopeTreeL2="",
             SCopeTreeL3="")

lp.create(filename="test3_endo_pyscenic_integrated_output.loom",
          layers=lf[:, :],
          row_attrs=row_attrs,
          col_attrs=col_attrs,
          file_attrs=attrs)

# Clustering the data
import seaborn as sns

# standard_scale ==> 0: standardize rows (0) or columns (1). By "standardize", it means for each row or column, subtract the minimum and divide each by its maximum
sns.clustermap(auc_mtx, method='ward', metric='euclidean', z_score=None, cmap="rocket_r")
sns.clustermap(auc_mtx, method='ward', metric='euclidean', z_score=None, standard_scale=1, cmap="rocket_r",
               cbar_pos=(0, 0, 0, 0))

#
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import loompy as lp
import json
import zlib
import scanpy.external as sce
from pyscenic.rss import regulon_specificity_scores
from pyscenic.cli.utils import load_signatures
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
import pandas as pd

test3_endo = sc.read_h5ad("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/test3_endo.h5ad")
test3_endo = test3_endo[test3_endo.obs['endo_leiden_r05'].isin(['0', '1', '2', '3'])]

sig = load_signatures(
    '/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/pySCENIC/EC_new2/test3_endo_reg.csv')
test3_endo = add_scenic_metadata(test3_endo, auc_mtx, sig)  # AUCell score가 test3_endo에 추가된다.

# 이 아래는 AUCell score가 추가되었을 때 X_aucell을 가지고 umap projection을 하는 방법.
# sce.pp.bbknn(test3_endo, batch_key='batch', use_rep='X_aucell')
# sc.tl.umap(test3_endo)
# sc.pl.umap(test3_endo, color=['batch'], add_outline=False, legend_loc='right margin', size=150, color_map='CMRmap')


from pyscenic.rss import regulon_specificity_scores
endo_leiden_to_celltype_dict = {'0': 'EC_1',
'1': 'EC_4',
'2': 'EC_2',
'3': 'EC_3'}

test3_endo.obs['EC_subclusters'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endo_leiden_to_celltype_dict[x]).astype('category')

reordered = ('EC_1', 'EC_2', 'EC_3', 'EC_4')
test3_endo.obs['EC_subclusters'] = test3_endo.obs['EC_subclusters'].cat.reorder_categories(list(reordered), ordered=True)

cellAnnot = test3_endo.obs[['batch', 'EC_subclusters']]
rss_EC_subclusters = regulon_specificity_scores(auc_mtx, cellAnnot['EC_subclusters'])


fig = plt.figure(figsize=(15, 8))

for c,num in zip(cats, range(1,len(cats)+1)):
    x = rss_EC_subclusters.T[c]
    ax = fig.add_subplot(1,4,num)
    plot_rss(rss_EC_subclusters, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )


rss_EC_subclusters.loc['EC_1'].sort_values(ascending=False).index[:20]
test3_endo.var[test3_endo.var['Regulon(Tbx3(+))']].index

# 아래처럼 Kdm4a regulon은 '1'에 높이 나타나는데 (RSS: 0.5260626882997943) 실제 plot에 expression을 놓고 보면 그렇지 않아 보임. 이것은 test3_endo의 raw counts를 가지고 했기 때문일까?
# 즉, size_factor normalization을 하고 난 뒤의 count를 pySCENIC input으로 집어넣어야할까?
test3_endo.var[test3_endo.var['Regulon(Kdm4a(+))']].index
sc.pl.umap(test3_endo, color=['Kdm4a', 'Tgfb2'], layer='magic', color_map=cmap)

# CELL TYPE SPECIFIC REGULATORS - Z-SCORE

signature_column_names = list(test3_endo.obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda x: x.startswith('Regulon('), signature_column_names))
test3_endo_scores = test3_endo.obs[signature_column_names + ['EC_subclusters']]
test3_endo_results = ((test3_endo_scores.groupby('EC_subclusters').mean() - test3_endo.obs[signature_column_names].mean()) / test3_endo.obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0: 'Z'})
test3_endo_results['regulon'] = list(map(lambda x: x[8:-1], test3_endo_results.regulon))

signature_column_names = list(test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC_5', 'EC_6'])].obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda x: x.startswith('Regulon('), signature_column_names))
test3_endo_scores = test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC_5', 'EC_6'])].obs[signature_column_names + ['EC_subclusters']]
test3_endo_results = ((test3_endo_scores.groupby('EC_subclusters').mean() - test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC_5', 'EC_6'])].obs[signature_column_names].mean()) / test3_endo[~test3_endo.obs['EC_subclusters'].isin(['EC_5', 'EC_6'])].obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0: 'Z'})

# 위를 Z-score가 아니라 일반 scale로 하는 것도 방법

test3_endo_results['regulon'] = list(map(lambda x: x[8:-1], test3_endo_results.regulon))


test3_endo_heatmap = pd.pivot_table(data=test3_endo_results[(test3_endo_results.Z >= 1.25)].sort_values('Z', ascending=False), index='EC_subclusters', columns='regulon', values='Z')

ec1_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_1'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)
ec2_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_2'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)
ec3_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_3'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)
ec4_regulon_list = list(test3_endo_results[test3_endo_results['EC_subclusters'] == 'EC_4'].sort_values('Z', ascending=False).iloc[:10]['regulon'].values)

test3_endo_heatmap = pd.pivot_table(data=test3_endo_results[test3_endo_results['regulon'].isin(ec1_regulon_list + ec2_regulon_list + ec3_regulon_list + ec4_regulon_list)].sort_values('Z', ascending=False), index='EC_subclusters', columns='regulon', values='Z')

import matplotlib
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

sns.clustermap(test3_endo_heatmap.sort_index(ascending=False),
               row_cluster=False,
               method='ward',
               metric='euclidean', z_score=None, standard_scale=None, cmap=cmap, xticklabels=True, yticklabels=True)

