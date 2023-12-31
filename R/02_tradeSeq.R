library(tradeSeq)
BiocParallel::register(BiocParallel::SerialParam())

sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)
plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)



data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")
set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
									 nGenes = 200, verbose = T)
