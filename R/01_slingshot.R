library(Seurat)
library(cowplot)
library(slingshot)
infile="tmp/tmp_seurat.rds"
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

data <- readRDS(infile)
# Save the objects as separate matrices for input in slingshot
dimred <- data@reductions$umap@cell.embeddings
clustering <- data$RNA_snn_res.0.8
counts <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
curves <- getCurves(lineages, 
										approx_points = 300, 
										thresh = 0.01, 
										stretch = 0.8, 
										allow.breaks = FALSE, 
										shrink = 0.99)

sce <- fitGAM(counts = as.matrix(counts), sds = curves)
plotGeneCount(curves, counts, clusters = clustering, models = sce)

