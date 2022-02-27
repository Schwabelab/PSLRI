# 2020 Sep 09 created 
# this script of PSLRI integrates Proteomics (LCMS) data to scRNA-seq for Ligand-receptor interaction using CellPhoneDB framework


setwd("D:/pCloud Sync/Columbia/SchwabeLab/P2Y14/onlineShared/PSLRI")

library(Seurat)
library(readxl)

#************************************************
#************************************************ get the data
#************************************************
mouseWLiverNormal_RS025 <- readRDS(file = "data/seurat4_musWLiv_RS025.rds")

DimPlot(mouseWLiverNormal_RS025, label = T)
DimPlot(mouseWLiverNormal_RS025,group.by = "cellTypesRes4", label = T)
DimPlot(mouseWLiverNormal_RS025,group.by = "orig.ident", label = T)


dim(mouseWLiverNormal_RS025@assays$RNA@counts)
dim(mouseWLiverNormal_RS025@assays$RNA@data)
# get the raw gene expression counts
rawCounts <- as.matrix(mouseWLiverNormal_RS025@assays$RNA@counts)
# make the cell names in raw counts and normalized counts to match to map the cell type annotation
temp_pos <- match(colnames(mouseWLiverNormal_RS025@assays$RNA@data),colnames(mouseWLiverNormal_RS025@assays$RNA@counts))
colnames(mouseWLiverNormal_RS025@assays$RNA@data)[1:5]
colnames(mouseWLiverNormal_RS025@assays$RNA@counts)[temp_pos[1:5]]
rawCounts <- rawCounts[,temp_pos]
rawCounts[1:5,1:5]
# get the cell type annotation
cellLabels <- data.frame(cells=colnames(mouseWLiverNormal_RS025@assays$RNA@data), cellType=mouseWLiverNormal_RS025@meta.data$cellTypesRes4)
head(cellLabels)
dim(cellLabels)
all(colnames(rawCounts)==as.character(cellLabels[,"cells"]))
table(cellLabels[,"cellType"])


#reading in the metabolies, their logFC, and interactions
metabolites <- as.data.frame(read_excel("MetaboliteLigandReceptorList.xlsx", sheet = "metabolites", col_names = T))
metabolites

interactions <- as.data.frame(read_excel("MetaboliteLigandReceptorList.xlsx", sheet = "interaction", col_names = T))[,c(1,2)]
interactions
# 
#************************ integrating proteomics and scRNA-seq data
# create new cell types as dying cells which are assumed to be 10% of the total number of hepatocytes
# add the proteiomic measurement of the metabolites to the dying cells along with other expressed genes in scRNA-seq

# TO DO
# # Expression distribution of raw counts of top 50 genes in all cells
# temp_rawCounts <- rawCounts
# total.reads <- rowSums(temp_rawCounts)#
# sorted <- order(total.reads, decreasing = TRUE )#based on total reads
# temp_rawCounts[sorted[1:10],1:10]
# boxplot(t(temp_rawCounts[sorted[1:50],]), horizontal = TRUE, las =1,cex.axis=0.6, cex.lab=0.1)
# title(paste0("Expression distribution of raw counts of top 50 genes in all cells"))
# quantile(t(temp_rawCounts[sorted[1:50],]))
# # Expression distribution of raw counts of top 50 genes in Hepatocytes
# temp_pos <- which(cellLabels[,"cellType"]=="Hepatocytes")
# # temp_pos <- temp_pos[1:(0.1*length(temp_pos))] #get 10% of the cells
# temp_rawCounts <- rawCounts[,temp_pos] # temp_rawCounts <- rawCounts
# total.reads <- rowSums(temp_rawCounts)#
# sorted <- order(total.reads, decreasing = TRUE )#based on total reads
# temp_rawCounts[sorted[1:10],1:10]
# boxplot(t(temp_rawCounts[sorted[1:50],]), horizontal = TRUE, las =1,cex.axis=0.6, cex.lab=0.1)
# title(paste0("Expression distribution of raw counts of top 50 genes in Hepatocytes"))
# quantile(t(temp_rawCounts[sorted[1:50],]))

# checking the expression distribution of the receptors of the metabolities in all cells
temp_pos <- which(toupper(rownames(rawCounts))%in% interactions[,"Symbol"])
rawCounts[temp_pos,1:5]
quantile(as.matrix(rawCounts[temp_pos,]))
boxplot(t(as.matrix(rawCounts[temp_pos,])), horizontal = TRUE, las =1,cex.axis=0.6, cex.lab=0.1)
mean(as.matrix(rawCounts[temp_pos,]))
# 
#creating the metabolies as new gene-rows in the existing single-cell data
temp_mat <- matrix(data = 0, nrow = nrow(metabolites), ncol = ncol(rawCounts))
rownames(temp_mat) <- metabolites[,"New_name"]
colnames(temp_mat) <- colnames(rawCounts)
temp_mat[,1:5]
dim(temp_mat)
# temp_mat[,temp_pos] <- 4 #giving the raw counts of metabolites
# temp_mat[,1:5]
# temp_mat[,temp_pos[1:5]]

# adding the new rows
rawCounts_metab <- rbind(rawCounts,temp_mat)
dim(rawCounts)
dim(rawCounts_metab)
rawCounts_metab[1:5,1:5]
rawCounts_metab[(nrow(rawCounts_metab)-5):nrow(rawCounts_metab),1:5]
# rawCounts_metab[(nrow(rawCounts_metab)-5):nrow(rawCounts_metab),temp_pos[1:5]]

#creating the new dying cells' columns
temp_pos <- which(cellLabels[,"cellType"]=="Hepatocytes")
temp_pos <- temp_pos[1:(0.1*length(temp_pos))] #get 10% of the cells
dyingCells <- rawCounts[,temp_pos]
dim(dyingCells)
colnames(dyingCells) <- paste0(colnames(dyingCells), "_dCells")
dyingCells[1:10,1:10]
temp_mat <- matrix(data = 0, nrow = nrow(metabolites), ncol = ncol(dyingCells))
rownames(temp_mat) <- metabolites[,"New_name"]
colnames(temp_mat) <- colnames(dyingCells)
temp_mat[,1:5]
# quantile(as.matrix(dyingCells))
# mean(as.matrix(dyingCells))
# quantile(as.matrix(rawCounts))
# mean(as.matrix(rawCounts))
temp_mat[,] <- max(1,median(as.matrix(rawCounts))) # giving the initial raw counts of metabolites
temp_mat[,1:5]
temp_mat <- temp_mat * metabolites[,"Fold_change"] # scaling the initial raw counts by fold-change to capture their expression variations
temp_mat[,1:5]
temp_mat <- round(temp_mat)
temp_mat[,1:5]
dyingCells_metab <- rbind(dyingCells,temp_mat)
dim(dyingCells)
dim(dyingCells_metab)
dyingCells_metab[1:5,1:5]
dyingCells_metab[(nrow(dyingCells_metab)-5):nrow(dyingCells_metab),1:5]
dyingCells_metab[(nrow(dyingCells_metab)-5):nrow(dyingCells_metab),1:5]

#adding the dying cells to other cells in single-cell data
rawCounts_metab <- cbind(rawCounts_metab,dyingCells_metab)
dim(rawCounts)
dim(rawCounts_metab)
rawCounts_metab[1:5,1:5]
rawCounts_metab[(nrow(rawCounts_metab)-5):nrow(rawCounts_metab),1:5]
rawCounts_metab[(nrow(rawCounts_metab)-5):nrow(rawCounts_metab),(ncol(rawCounts_metab)-5):ncol(rawCounts_metab)]

# normalizing for library depth
temp <- rawCounts_metab #as.matrix(mouseLiverCCA_seurat2_rs030@assays$RNA@counts)
dim(temp)
temp[1:5,1:5]
temp_norm_mat <- apply(temp, 2, function(x) (x/sum(x))*10000) # cpm normalization to 1e4
colSums(temp_norm_mat[,1:10])
rowSums(temp_norm_mat)[c(1:10,(nrow(rawCounts_metab)-5):nrow(rawCounts_metab))]
temp_norm_mat[(nrow(rawCounts_metab)-nrow(metabolites)+1):nrow(rawCounts_metab),1:5]

#getting the metabolites separately to add after mouse to human gene conversion since CellphneDB works only with human genes
temp_norm_mat_metab <- temp_norm_mat[(nrow(rawCounts_metab)-nrow(metabolites)+1):nrow(rawCounts_metab),]
temp_norm_mat_metab[,1:5]
temp_norm_mat_metab[,(ncol(temp_norm_mat_metab)-5):ncol(temp_norm_mat_metab)]
dim(temp_norm_mat_metab)
remove(temp)
# 
# 
#doing the mouse to human gene converstion by getting the homologs from MGI JAX
humGenes_req <- readRDS(file = "data/homologs_mouse_human_geneSymbol.rds")
dim(humGenes_req)
head(humGenes_req)


# converting the mouse GE matrix to human
#finding the matching genes in the matrix
dim(temp_norm_mat) #[1] 16089  2492
temp_norm_mat[1:5,1:5]
dim(humGenes_req)
temp_pos <- which(rownames(temp_norm_mat) %in% humGenes_req[,"MGI.symbol"])
# all(rownames(temp_norm_mat)[temp_pos]==humGenes_req[,"MGI.symbol"])
temp_norm_mat_req <- temp_norm_mat[temp_pos,]
dim(temp_norm_mat_req)
temp_norm_mat_req[1:5,1:5]
#get the required human genes present in mouse GE matrix
temp_pos <- which(humGenes_req[,"MGI.symbol"] %in% rownames(temp_norm_mat))
# all(rownames(temp_norm_mat)[temp_pos]==humGenes_req[,"MGI.symbol"])
humGenes_req_map <- humGenes_req[temp_pos,]
dim(humGenes_req_map)


# ordering the mouse and human genes appropriately for mapping
temp_pos <- match(rownames(temp_norm_mat_req),humGenes_req_map[,"MGI.symbol"])
humGenes_req_map[temp_pos,"MGI.symbol"][1:10]
rownames(temp_norm_mat_req)[1:10]
humGenes_req_map <- humGenes_req_map[temp_pos,]
all(rownames(temp_norm_mat_req)==humGenes_req_map[,"MGI.symbol"])
temp_norm_mat_req[1:5,1:5]
rownames(temp_norm_mat_req) <- humGenes_req_map[,"HGNC.symbol"]
temp_norm_mat_req[1:5,1:5]
dim(temp_norm_mat_req)
mus2humGenesNormData <- temp_norm_mat_req
dim(mus2humGenesNormData)
mus2humGenesNormData[1:5,1:5]

# adding back the metabolites
length(which(rownames(mus2humGenesNormData)%in%metabolites[,1])) #checking if the metabolites are conserved (it should not be)
mus2humGenesNormData_metab <- rbind(mus2humGenesNormData, temp_norm_mat_metab)
dim(mus2humGenesNormData)
dim(mus2humGenesNormData_metab)
mus2humGenesNormData_metab[1:5,1:5]
mus2humGenesNormData_metab[(nrow(mus2humGenesNormData_metab)-nrow(metabolites)+1):nrow(mus2humGenesNormData_metab),1:5]
mus2humGenesNormData_metab[(nrow(mus2humGenesNormData_metab)-nrow(metabolites)+1):nrow(mus2humGenesNormData_metab),(ncol(mus2humGenesNormData_metab)-5):ncol(mus2humGenesNormData_metab)]


# #**********************************adding the dying cell annotation to the current cell type annotation
#add the dying cells
head(cellLabels)
dim(cellLabels)
dim(mus2humGenesNormData_metab)
temp_pos <- match(cellLabels[,"cells"], colnames(mus2humGenesNormData_metab))
temp <- data.frame(cells=colnames(mus2humGenesNormData_metab)[-temp_pos],cellType=rep("Dying cell", length(colnames(mus2humGenesNormData_metab)[-temp_pos])))
head(temp)
cellLabels_req <- rbind(cellLabels, temp)
dim(cellLabels_req)
dim(mus2humGenesNormData_metab)
all(colnames(mus2humGenesNormData_metab)==cellLabels_req[,"cells"])

# removing multi clusters and blood cells
table(cellLabels_req[,"cellType"])
#getting the required cell types
temp_pos <- which(cellLabels_req[,"cellType"]%in% c("Erythroid", "Multi"))
cellLabels_req <- cellLabels_req[-temp_pos, ]
head(cellLabels_req)
table(cellLabels_req[,"cellType"])

# getting the required cell counts
temp_pos <- match(cellLabels_req[,"cells"], colnames(mus2humGenesNormData_metab))
cellLabels_req[1:5,"cells"]
colnames(mus2humGenesNormData_metab)[temp_pos[1:5]]
mus2humGenesNormData_metab_req <- mus2humGenesNormData_metab[,temp_pos]
mus2humGenesNormData_metab_req[1:5,1:5]
all(cellLabels_req[,"cells"]==colnames(mus2humGenesNormData_metab_req))
dim(mus2humGenesNormData_metab_req)
dim(cellLabels_req)

#**********************************getting the expression data for CellphoneDB
write.table(t(c("Genes", colnames(mus2humGenesNormData_metab_req))),file = "results/data_cellphoneDB_countNorm_musRS0025_mgiJax_metab_4_newMetabATP.txt", quote = FALSE, sep = "\t", append = FALSE, col.names = FALSE, row.names = FALSE)
write.table(mus2humGenesNormData_metab_req,file = "results/data_cellphoneDB_countNorm_musRS0025_mgiJax_metab_4_newMetabATP.txt", quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE)


# #getting the metadata
head(cellLabels_req)
table(cellLabels_req$cellType)
dim(cellLabels_req)
dim(mus2humGenesNormData_metab_req)

metaData <- data.frame(Cell=cellLabels_req[,"cells"],cell_type=cellLabels_req[,"cellType"])
head(metaData)
metaData$cell_type <- as.character(metaData$cell_type)
table(metaData[,"cell_type"])
dim(metaData)

# changing B and plasma cells
# metaData$cell_type <- as.character(metaData$cell_type)
tmp_pos <- which(metaData$cell_type=="B cell")
metaData[tmp_pos,]
metaData$cell_type[tmp_pos] <- "Plasma cell"
metaData[tmp_pos,]
tmp_pos <- which(metaData$cell_type=="Liver resident B cell")
metaData[tmp_pos,]
metaData$cell_type[tmp_pos] <- "B cell"
metaData[tmp_pos,]
table(metaData[,"cell_type"])
head(metaData)

all(metaData[,"Cell"]==colnames(mus2humGenesNormData_metab_req))

#**********************************getting the metadata for CellphoneDB
write.table(metaData, file = "results/data_cellphoneDB_metaData_musRS0025_4_newMetabATP.txt", sep = '\t', quote = F, row.names = F, col.names = T)




