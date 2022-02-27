# 2020 Sep 10
# the script to run CellPhoneDB with the new database


# Activate virtual-env
source cpdb-venv/bin/activate

# go to the program folder
cd /Users/an2863/pCloud\ Sync/Columbia/SchwabeLab/P2Y14/cpdb/

# run CellPhoneDB
cellphonedb method statistical_analysis data_cellphoneDB_metaData_musRS0025_4_newMetabATP.txt data_cellphoneDB_countNorm_musRS0025_mgiJax_metab_4_newMetabATP.txt --output-path=out_newDBUser_metab_4_newMetab_ATP_RS025 --counts-data=gene_name --database out/cellphonedb_user_2021-12-16-15_41.db