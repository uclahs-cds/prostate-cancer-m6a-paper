library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);


sample.clusters <- read.delim(
	file = '2022-01-18_m6a_adjusted_peak_cluster_assignment.txt',
	as.is = TRUE
	);

clusters.ordered <- sample.clusters[order(sample.clusters$cluster), ];

muts <- read.delim('2020-02-25_CPC-GENE_clinical_from_database.txt',
	as.is = TRUE
	);

patient.overlap <- intersect(clusters.ordered$sample, muts$sample_id);

m6a.clusters <- clusters.ordered[match(patient.overlap, clusters.ordered$sample), ];
clinical <- muts[match(patient.overlap, muts$sample_id), c('sample_id', 'pga')];

forplot <- data.frame(
	cluster = as.factor(m6a.clusters$cluster),
	pga = clinical$pga
	);

aov.results <- aov(
	formula = pga ~ cluster,
	data = forplot
	);

pvalue <- summary(aov.results)[[1]][['Pr(>F)']][1];

your.number <- list(
	base = scientific.notation(pvalue, digits = 2, type = 'list')$base,
	exponent = scientific.notation(pvalue, type = 'list')$exponent
	);

stats.key <- list(
	text = list(
		lab =  as.expression(substitute(P == paste(base %*% 10^exponent),  your.number)),
		cex = .6,
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
		x = 0.01,
		y = 0.98,
		corner = c(0, 1)
		)
	);

col <- c('royalblue2', 'firebrick', 'chartreuse3', 'purple4', 'darkorange')[match(forplot$cluster, c(1, 2, 3, 4, 5))];
create.boxplot(
	formula = pga ~ cluster,
	height = 2.5,
	width = 2.75,
	data = forplot,
	filename = 'SupplementaryFigure3A.png',
	resolution = 300,
	add.stripplot = TRUE,
	ylab.label = 'PGA',
	ylab.cex = .6,
	xlab.cex = 0,
	points.col = col,
	points.cex = 0.25,
	xaxis.lab = c('P1', 'P2', 'P3', 'P4', 'P5'),
	legend = stats.legend,
	xaxis.tck = c(.6, 0),
	yaxis.tck = c(.6, 0),
	xaxis.cex = .6,
	yaxis.cex = .6,
	);
