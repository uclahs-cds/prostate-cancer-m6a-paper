library(BoutrosLab.plotting.general)
library(cluster)
library(apcluster)

# dissimilarity matrix is 1 - simmilarity matrix
load(
	file = '2021-09-28_adjusted_peak_sample_cluster_results.rda'
	)

sample.cluster.n <- 5

class <- sample.cluster[[sample.cluster.n]]$consensusClass
dmatrix <- 1 - sample.cluster[[sample.cluster.n]]$consensusMatrix
rownames(dmatrix) <- colnames(dmatrix) <- names(class)
score <- silhouette(
	x = class,
	dmatrix = dmatrix
	)
pdf('SupplementaryFigure5Ci.pdf')
plot(score)
dev.off()

# and gene clusters
load(
	file = '2021-10-05_adjusted_peak_gene_cluster_results.rda'
	)

gene.cluster.n <- 5

gene.class <- gene.cluster[[gene.cluster.n]]$consensusClass
gene.dmatrix <- 1 - gene.cluster[[gene.cluster.n]]$consensusMatrix
gene.score = silhouette(gene.class, dmatrix = gene.dmatrix)

pdf('SupplementaryFigure5Cii.pdf')
plot(gene.score)
dev.off()
