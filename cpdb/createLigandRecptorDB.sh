# 2020 Sep 10
# the script to create a new Ligand-receptor interaction pair database in CellPhoneDB

# Activate CellPhoneDB virtual-env
source cpdb-venv/bin/activate

# go to the program folder
cd /Users/an2863/pCloud\ Sync/Columbia/SchwabeLab/P2Y14/cpdb/

#generate new DB
cellphonedb database generate --user-interactions newInteractions.csv --user-protein newProteins.csv --user-gene newGenes.csv