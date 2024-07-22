library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

sample.clusters <- read.delim(
	file = '2022-01-18_m6a_adjusted_peak_cluster_assignment.txt',
	as.is = TRUE
	);

clusters.ordered <- sample.clusters[order(sample.clusters$cluster), ];

muts <- read.delim(
	file = '2020-02-25_CPC-GENE_clinical_from_database.txt',
	as.is = TRUE
	);

patient.overlap <- intersect(clusters.ordered$sample, muts$sample_id);

m6a.clusters <- clusters.ordered[match(patient.overlap, clusters.ordered$sample), ];
clinical <- muts[match(patient.overlap, muts$sample_id), c('sample_id', 'CellularityPath')];

table(clinical$CellularityPath, m6a.clusters$cluster)

forplot <- data.frame(
	x = as.factor(m6a.clusters$cluster),
	y = clinical$CellularityPath
	);

aov.results <- aov(
	formula = y ~ x,
	data = forplot
	);

pvalue <- summary(aov.results)[[1]][['Pr(>F)']][1];

stats.key <- list(
	text = list(
		lab = paste('P =', round(pvalue, digits = 3)),
		cex = 1.25,
		fontface = 'bold'
		),
	padding.text = 0,
	x = 0.06,
	y = 0.97,
	corner = c(0, 1),
	title = NULL,
	cex.title = 1,
	background = 'white',
	alpha.background = 0
	);

stats.legend <- list(
	inside = list(
		fun = draw.key,
		args = list(
			key = stats.key
			),
		x = 0.75,
		y = 0.98
		)
	);

create.boxplot(
	formula = y ~ x,
	data = forplot,
	filename = generate.filename('m6A', 'SupplementaryFigure5E', 'png'),
	resolution = 300,
	ylab.label = c('Tumour Purity'),
	ylab.cex = 1,
	xlab.cex = 0,
	xaxis.lab = c('P1', 'P2', 'P3', 'P4', 'P5'),
	xaxis.cex = 1,
	yaxis.cex = 1,
	);
