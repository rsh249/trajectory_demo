library(Seurat)
library(cowplot)
library(slingshot)

if(!dir.exists("tmp")){dir.create("tmp")}


args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
if(is.na(infile)){infile="tmp/tmp_seurat.rds"}



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

pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

## save and run fitGAM as it's own process for parallelization
write.table(pseudotime, file="tmp/tmp_pseudo.mat", sep="/t")
write.table(cellWeights, file="tmp/tmp_weights.mat", sep="/t")
saveRDS(curves, file="tmp/tmp_curves.rds")
saveRDS(lineages, file="tmp/tmp_lineages.rds")


