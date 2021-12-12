# Get mouse ortholog of human genes using biomaRt
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

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


cellphonedb plot dot_plot \
--means-path means.txt \
--pvalues-path pvalues.txt \
--output-path ./Figure \
--output-name dotplot.pdf \
--verbose

cellphonedb plot heatmap_plot \
../../test_meta.txt \
--pvalues-path pvalues.txt \
--output-path ./Figure \
--count-name heatmap_count.pdf \
--log-name heatmap_log_count.pdf \
--count-network-name test_count_network.txt \
--interaction-count-name test_interactions_count.txt

