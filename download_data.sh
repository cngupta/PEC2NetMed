# Download GRNs and other data to reproduce analysis as shown in paper

mkdir -p data/
mkdir -p tables/ # for outputs
mkdir -p data/GNN_labels
mkdir -p data/gtrddb

# links to drugbank files are not provided here.
# Please acquire a licence by writing to them.
# Two files required in the data/ to reproduce the analysis in PEC2 NetMed paper.
# 1) drugbank.tsv
# 2) pharmacologically_active.csv

# links to disgenet db are not provided.
# Please acquire a licence at https://www.disgenet.org/downloads
# The all_gene_disease_associations.tsv file is required in  data/


# DGIB data
wget -P https://www.dgidb.org/data/latest/interactions.tsv

#GRNs
wget -P data/ "https://zenodo.org/records/14252962/files/PEC2_NetMed_GRNs.zip"
cd data/
unzip PEC2_NetMed_GRNs.zip

# GTRD database
cd gtrddb
# create a files.txt with list of all TF IDs you want. Check GTRD db website
while read file; do
  curl -O "http://gtrd.biouml.org:8888/downloads/current/intervals/target_genes/Homo%20sapiens/genes%20whole%5b-5000,+5000%5d/$file"
done < files.txt
cd ../../

# PEC2 DEGs
wget -P data https://brainscope.gersteinlab.org/data/DEG-combined/ASD_DEGcombined.csv
wget -P data https://brainscope.gersteinlab.org/data/DEG-combined/Bipolar_DEGcombined.csv
wget -P data https://brainscope.gersteinlab.org/data/DEG-combined/Schizophrenia_DEGcombined.csv

# GO genesets
wget -P data https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/GO/Human_GO_bp_with_GO_iea_symbol.gmt

# LINCS
wget -P data https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gmt/l1000_cp.gmt
wget -P data https://s3.amazonaws.com/lincs-dcic/sigcom-lincs-metadata/LINCS_small_molecules.tsv
wget -P data https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/means/cp_mean_coeff_mat.tsv.gz
cd data
gunzip cp_mean_coeff_mat.tsv.gz
cd ../

# B2DB
wget -P data/ "https://figshare.com/ndownloader/files/30561027"
cd data
unzip 30561027 # this might change based on date; its a zipped folder
mv B3DB-main/B3DB/B3DB_classification.tsv .
rm -r 30561027 B3DB-main
cd ../
