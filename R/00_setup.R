library(Seurat)
library(glue)

args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
if(is.na(infile)){infile="data/GSE72857_umitab.txt"} ## keep for manual testing of this R script
if(!dir.exists("tmp")){dir.create("tmp")}
outSeurat="tmp/tmp_seurat.rds"

message(glue("Converting gene counts {infile} to matrix"))
data <- read.delim(infile, header = T, row.names = 1)
comp_matrix <- Matrix::Matrix(as.matrix(data), sparse = T)
message(glue("Running basic Seurat pipeline for dimensionality reduction and clustering"))


umi_counts <- comp_matrix
# Data analysis with Seurat pipeline
data <- CreateSeuratObject(counts = umi_counts)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data,resolution = 1)
data <- FindClusters(data)
data <- RunUMAP(data, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)

message(glue("Saving Seurat object as {outSeurat}"))
saveRDS(data, outSeurat)
