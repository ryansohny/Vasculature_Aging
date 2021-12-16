# Get mouse ortholog of human genes using biomaRt
# Get mouse ortholog of human genes using biomaRt
library(biomaRt)

mart <- useMart('ENSEMBL_MART_ENSEMBL')
# > listMarts()
#                biomart                version
# 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 105
# 2   ENSEMBL_MART_MOUSE      Mouse strains 105
# 3     ENSEMBL_MART_SNP  Ensembl Variation 105
# 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 105
mart <- useDataset('hsapiens_gene_ensembl', mart)
# listDatasets(mart) ==> mart에 존재하는 여러 datasets 중에서 'hsapiens_gene_ensembl을 쓰는 것

annotLookup <- getBM(
    mart = mart,
    attributes = c('ensembl_gene_id',
                   'gene_biotype',
                   'external_gene_name',
                   'uniprot_gn_symbol',
                   'uniprot_gn_id',
                   'hgnc_symbol'),
    uniqueRows=TRUE)

# source activate cellphonedb
# /home/phenomata/cpdb/releases/v2.0.0/data
# /data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB

#####################################################################################################################################
#####################################################################################################################################

# From https://www.ensembl.org/biomart/martview/e76ba4ad9601a50ce208f4f2a1742513
# Dataset: Ensembl Genes 104 ==> Mouse genes (GRCm39)
# ==> Attributes: Gene stable ID, Gene stable ID version, Transcript stable ID, Gene name, Human gene stable ID, Human gene name and Human homology type
# tsv download and pick only "ortholog_one2one" in Human homology type field
# ==> /data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB/Ensembl_biomaRt_hs-mm_ortholog_one2one.txt

#####################################################################################################################################
#####################################################################################################################################


## only Ortholog gene_input generation
dbf = open("/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB/Ensembl_biomaRt_hs-mm_ortholog_one2one.txt", 'r')
db = list(x.strip('\n').split('\t')[3] for x in dbf)
db.remove('Human gene stable ID')
dbf.close()

dfh = open("gene_input.csv", 'r')
rfh = open("gene_input_hs-mm_ortholog.csv", 'w')
rfh.write(dfh.readline())
for i in dfh:
    if i.strip('\n').split(',')[3] in db:
        rfh.write(i)

## Custom Database Generation (i)
cellphonedb database generate \
--user-protein protein_curated.csv \
--user-gene gene_input_hs-mm_ortholog.csv \
--user-complex complex_curated.csv \
--user-interactions interaction_curated.csv \
--user-interactions-only \
--result-path ./CustomDatabase_20211216_1 \
--log-file customdatabase_20211216_1.out

## Custom Database Generation (ii)
cellphonedb database generate \
--user-protein protein_curated.csv \
--user-gene gene_input_hs-mm_ortholog.csv \
--user-complex complex_curated.csv \
--user-interactions interaction_curated.csv \
--user-interactions-only \
--fetch \
--result-path ./CustomDatabase_20211216_2 \
--log-file customdatabase_20211216_2.out

#####################################################################################################################################
#####################################################################################################################################

cellphonedb method statistical_analysis \
--project-name test3_stat \
--output-path test3_stat \
--counts-data gene_name \
--threads 20 \
test3_meta.txt \
test3_forCellPhoneDB_uppercase.txt

cellphonedb method statistical_analysis \
--project-name test3_stat2 \
--output-path test3_stat2 \
--counts-data gene_name \
--threads 20 \
test3_meta2.txt \
test3_forCellPhoneDB_uppercase.txt

# threshold 10%
nohup cellphonedb method statistical_analysis \
--project-name test3_stat3 \
--output-path test3_stat3 \
--counts-data gene_name \
--threads 20 \
--threshold 0.1 \
test3_meta.txt \
test3_forCellPhoneDB_uppercase.txt > stat3.out &

# threshold 10% and pvalue < 0.005
nohup cellphonedb method statistical_analysis \
--project-name test3_stat4 \
--output-path test3_stat4 \
--counts-data gene_name \
--threads 20 \
--threshold 0.1 \
--pvalue 0.005 \
test3_meta.txt \
test3_forCellPhoneDB_uppercase.txt > stat4.out &


# threshold 10% and pvalue < 0.005 (Using Custom database)
/data/Projects/phenomata/01.Projects/11.Vascular_Aging/Database/CellPhoneDB/CustomDatabase_20211216/cellphonedb_user_2021-12-16-17_07.db












# only using means.txt
cellphonedb plot dot_plot \
--means-path means.txt \
--pvalues-path pvalues.txt \
--output-path ./Figure \
--output-name dotplot.pdf \
--verbose

cellphonedb plot heatmap_plot \
../../test3_meta.txt \
--pvalues-path pvalues.txt \
--output-path ./Figure \
--count-name heatmap_count.pdf \
--log-name heatmap_log_count.pdf \
--count-network-name test_count_network.txt \
--interaction-count-name test_interactions_count.txt
