# PSLRI
PSLRI (Proteomics scRNA-seq ligand receptor interaction) is a method to (a) integrate the expression of ligands and their receptors from proteomics (LCMS) and scRNA-seq data, and (b) to infer the ligand-receptor interaction between cell-types in the combined dataset. In the example provided here, the ligands expressed by dying cells are measured by proteomics. The gene expression of the different cell-types in a diseased mouse liver (CCl4 treatment) is measured using scRNA-seq. The proteomics measurement of the dying cell ligands are integrated to scRNA-seq data of dying Hepatocytes after appropriate normalization and scaling. After the proteomics and scRNA-seq are integrated, the ligand-receptor interactions are inferred using the CellPhoneDB framework. For this, the interacting ligand-receptor pairs are first added to the CellPhoneDB database because they were not already present. Then, ligand-receptor interaction inference is performed using CellPhoneDB.

The latest updates may be found at the original GitHub page: https://github.com/Schwabelab/PSLRI

Large input Seurat object file of the mouse whole-liver is at figshare page: https://figshare.com/articles/dataset/PSLRI_Proteomic_scRNA-seq_ligand-receptor_interaction_/19245966

DOI: https://doi.org/10.6084/m9.figshare.19245966.v2


<a href="https://doi.org/10.5281/zenodo.6301768"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6301768.svg" alt="DOI"></a>
