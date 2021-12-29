# Distinguishing markers for [0,2] – [1,3]
matplotlib.use('TkAgg')
%matplotlib

endoleiden_to_new_dict = {'0':'Vcam1-high', '1':'Cd36-high', '2':'Vcam1-high', '3':'Cd36-high', '4':'NA', '5':'NA'}

test3_endo.obs['twogroups'] = test3_endo.obs['endo_leiden_r05'].map(lambda x: endoleiden_to_new_dict[x]).astype('category')

sc.pl.umap(test3_endo, color=['twogroups', 'batch'], add_outline=False, legend_loc='right margin', size=130, color_map='CMRmap')

# Wilcoxon rank-sum test (Vcam1-high > Cd36-high)
sc.tl.rank_genes_groups(test3_endo, 'twogroups', method='wilcoxon', pts=True, key_added='Vcam1-highvsCd36-high_rank_genes_groups', groups=['Vcam1-high'], reference='Cd36-high')
sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['twogroups'].isin(['Vcam1-high', 'Cd36-high'])], n_genes=40, groups=['Vcam1-high'], key='Vcam1-highvsCd36-high_rank_genes_groups', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis')
# Wilcoxon rank-sum test (Cd36-high > Vcam1-high)
sc.tl.rank_genes_groups(test3_endo, 'twogroups', method='wilcoxon', pts=True, key_added='Cd36-highvsVcam1-high_rank_genes_groups', groups=['Cd36-high'], reference='Vcam1-high')
sc.pl.rank_genes_groups_heatmap(test3_endo[test3_endo.obs['twogroups'].isin(['Vcam1-high', 'Cd36-high'])], n_genes=40, groups=['Cd36-high'], key='Cd36-highvsVcam1-high_rank_genes_groups', show_gene_labels=True, min_logfoldchange=2, dendrogram=False, cmap='cividis')

result = test3_endo.uns['Vcam1-highvsCd36-high_rank_genes_groups']
groups = result['names'].dtype.names
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})
deg_wilcoxon = deg_wilcoxon.set_index('Vcam1-high_names')

deg_wilcoxon = pd.concat([deg_wilcoxon, test3_endo.var['mean'], test3_endo.var['n_cells']], axis=1)
sb.histplot(deg_wilcoxon['Vcam1-high_logfoldchanges'], kde=True, stat='count')
sb.jointplot(data=deg_wilcoxon, x='Vcam1-high_logfoldchanges', y='mean')


#######################################################################################################################################
# Gene-set enrichment analysis using gseapy
# GO_Biological_Process_2021 (deg_wilcoxon['Vcam1-high_scores'] > 10)
# source activate scanpy_1.8.1
import gseapy as gp
gp.get_library_name(organism='Mouse') # list of available mouse gene-set
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_scores'] > 10) & (deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05)].index.tolist() # 287
enrich_GOBP_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['GO_Biological_Process_2021'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_GOBP_Vcam1_up_scorewise',
cutoff=0.05)

degs_up_Cd36_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_scores'] < -10) & (deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05)].index.tolist()
enrich_GOBP_Cd36_up = gp.enrichr(gene_list=degs_up_Cd36_high,
gene_sets=['GO_Biological_Process_2021'],
organism='Mouse',
description='degs_up_Cd36_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_GOBP_Cd36_up_scorewise',
cutoff=0.05)

# Gene-set enrichment analysis using gseapy
# GO_Biological_Process_2021
gp.get_library_name(organism='Mouse')
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] > 5)].index.tolist() # 287
enrich_GOBP_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['GO_Biological_Process_2021'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_GOBP_up',
cutoff=0.05)

degs_up_Cd36_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] < -5)].index.tolist()
enrich_GOBP_Cd36_up = gp.enrichr(gene_list=degs_up_Cd36_high,
gene_sets=['GO_Biological_Process_2021'],
organism='Mouse',
description='degs_up_Cd36_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_GOBP_Cd36_up',
cutoff=0.05)

# MSigDB_Hallmark_2020
gp.get_library_name(organism='Mouse')
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] > 5)].index.tolist() # 287
enrich_Hallmark_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['MSigDB_Hallmark_2020'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_Hallmark_Vcam1_up',
cutoff=0.05)

degs_up_Cd36_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] < -5)].index.tolist()
enrich_Hallmark_Cd36_up = gp.enrichr(gene_list=degs_up_Cd36_high,
gene_sets=['MSigDB_Hallmark_2020'],
organism='Mouse',
description='degs_up_Cd36_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_Hallmark_Cd36_up',
cutoff=0.05)

# 'KEGG_2021_Human' ==> 이거 mouse는 없나
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] > 5)].index.tolist() # 287
enrich_kegg_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['KEGG_2021_Human'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_kegg_Vcam1_up',
cutoff=0.05)

degs_up_Cd36_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] < -5)].index.tolist()
enrich_kegg_Cd36_up = gp.enrichr(gene_list=degs_up_Cd36_high,
gene_sets=['KEGG_2021_Human'],
organism='Mouse',
description='degs_up_Cd36_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_kegg_Cd36_up',
cutoff=0.05)

# 'KEGG_2019_Mouse'
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] > 5)].index.tolist() # 287
enrich_kegg_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['KEGG_2019_Mouse'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_kegg_Vcam1_up',
cutoff=0.05)

degs_up_Cd36_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] < -5)].index.tolist()
enrich_kegg_Cd36_up = gp.enrichr(gene_list=degs_up_Cd36_high,
gene_sets=['KEGG_2019_Mouse'],
organism='Mouse',
description='degs_up_Cd36_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_kegg_Cd36_up',
cutoff=0.05)

# 'Mouse_Gene_Atlas'
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] > 5)].index.tolist() # 287
enrich_mga_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['Mouse_Gene_Atlas'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_mga_Vcam1_up',
cutoff=0.05)

degs_up_Cd36_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] < -5)].index.tolist()
enrich_mga_Cd36_up = gp.enrichr(gene_list=degs_up_Cd36_high,
gene_sets=['Mouse_Gene_Atlas'],
organism='Mouse',
description='degs_up_Cd36_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_mga_Cd36_up',
cutoff=0.05)

#######################################################################################################################################
# Gene-set enrichment analysis using gseapy
result = test3_endo.uns['endo_leiden_r05_rank_genes_groups']
groups = result['names'].dtype.names
deg_wilcoxon = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'logfoldchanges', 'scores', 'pvals_adj']})
deg_wilcoxon = deg_wilcoxon.set_index('Vcam1-high_names')

deg_wilcoxon = pd.concat([deg_wilcoxon, test3_endo.var['mean'], test3_endo.var['n_cells']], axis=1)
sb.histplot(deg_wilcoxon['Vcam1-high_logfoldchanges'], kde=True, stat='count')
sb.jointplot(data=deg_wilcoxon, x='Vcam1-high_logfoldchanges', y='mean')

# GO_Biological_Process_2021
gp.get_library_name(organism='Mouse')
degs_up_Vcam1_high = deg_wilcoxon[(deg_wilcoxon['Vcam1-high_pvals_adj'] < 0.05) & (deg_wilcoxon['Vcam1-high_logfoldchanges'] > 5)].index.tolist() # 287
enrich_GOBP_Vcam1_up = gp.enrichr(gene_list=degs_up_Vcam1_high,
gene_sets=['GO_Biological_Process_2021'],
organism='Mouse',
description='degs_up_Vcam1_high',
outdir='/mnt/data/Projects/phenomata/01.Projects/11.Vascular_Aging/03.Scanpy/gseapy/enrich_DEGs_GOBP_up',
cutoff=0.05)




#######################################################################################################################################
# Transcription Factors
dfh = open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/pySCENIC/mm_mgi_tfs.txt", 'r')
mmtf = list(x.strip() for x in dfh.readlines())




a = list(test3_endo.obs['batch'].values)
b = list(test3_endo.obs['twogroups'].values)
c = list(map(lambda x: a[x] + '_' + b[x], list(range(len(a)))))
test3_endo.obs['aging_twogroups'] = c
sc.tl.rank_genes_groups(test3_endo, 'aging_twogroups', method='wilcoxon', pts=True, key_added='Vcam1-highvsCd36-high_rank_genes_groups_aging')





def scale_01(df):
tmp = df.values
row_index = df.index
col_index = df.columns
row_min = tmp.min(axis=1)
tmp =  tmp - row_min[:,None]
row_max = tmp.max(axis=1)
scaled = pd.DataFrame(tmp / row_max[:, None], index = row_index, columns = col_index)
return scaled


pval_df = test3_endo.obs.loc[~filter_cells,['ct_pseudotime']]
ind = np.argsort(pval_df.loc[:, 'ct_pseudotime'])

def GetGamTrends(pdt_cells, genes, data, n_splines = 8):
win = 30
min_val = 1e9
gamx, gamy, gam_ci = {},{},{}
for g in genes:
X = pdt_cells.values
y = np.array(data.loc[pdt_cells.index,g])
gam = LinearGAM(s(0, n_splines=n_splines,spline_order=3))
gam.fit(X, y)
gam.gridsearch(X, y, progress=False)
XX = gam.generate_X_grid(term=0)
m = X.min()
M = X.max()
XX = np.linspace(m - .10, M + .10, 500)
YY = gam.predict(XX)
gamx[g] = XX
gamy[g] = YY
if min_val > YY.min():
min_val = YY.min()
        gam_ci[g] = gam.prediction_intervals(XX, width=.95)
    return gamx, gamy, gam_ci


