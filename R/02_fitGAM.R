### FIT GAM
library(Seurat)
library(cowplot)
library(slingshot)
library(BiocParallel)


args <- commandArgs(trailingOnly = TRUE)
pseudo=args[1]
weights=args[2]
cvs=args[3]
lin=args[4]
nclus=args[5]
if(length(args)==0){
	pseudo="tmp/tmp_pseudo.mat"
	weights="tmp/tmp_weights.mat"
	cvs="tmp/tmp_curves.rds"
	lin="tmp/tmp_lineages.rds"
	nclus=2
}

pseudo=read.table(pseudo)
weights=read.table(weights)
curves=readRDS(cvs)
lineages=readRDS(lin)

## Parallelization
snowparam <- SnowParam(workers = nclus, type = "SOCK")
register(snowparam, default = TRUE)
registered()
##

sce <- fitGAM(counts = as.matrix(counts)[1:10,], 
							pseudotime = pseudotime, 
							cellWeights = cellWeights,
							nknots = 6, verbose = FALSE,
							parallel=T,
							BPPARAM=BiocParallel::bpparam())


plotGeneCount(curves, counts, clusters = clustering, models = sce)

