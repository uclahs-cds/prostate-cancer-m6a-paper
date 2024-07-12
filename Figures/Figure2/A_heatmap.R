library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(getopt);

source('heatmap_functions.R');

input <- 'adjusted'
level <- 'peak'
samplek <- 5
genek <- 5

meth <- read.delim(
	file = 'MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakCounts.Adjusted.tsv',
	as.is = TRUE
	);
load('2021-09-28_adjusted_peak_sample_cluster_results.rda');
load('2021-10-05_adjusted_peak_gene_cluster_results.rda');

# do the same filtering as before so you are working with the same matrix
binary <- ifelse(
	test = meth == 0,
	yes = 0,
	no = 1
	);
count <- rowSums(binary);
meth.filtered <- log2(meth[which(count == ncol(meth)), ]);

# also filter based on variance because there are too many peaks for consensus clustering
meth.var <- apply(
	X = meth.filtered,
	MARGIN = 1,
	FUN = IQR
	);
cutoff <- quantile(x = meth.var, probs = 0.75);
var.filtered <- meth.filtered[meth.var > cutoff, ];

# scale works on columns, so have to transpose to scale by genes
scale.meth <- t(scale(t(var.filtered)));

# read in any clinical information you have and want to add to the heatmap
clinical <- read.delim(
	file = '2020-02-25_CPC-GENE_clinical_from_database.txt',
	as.is = TRUE
	);
clinical$clinicalt3 <- tolower(substr(clinical$pathologic_t, 1, 2));

clinical.filtered <- clinical[
	match(colnames(meth), clinical$patient_id),
	c('age_at_treatment', 'pre_treatment_psa', 'clinicalt3', 'pathologic_isup_grade', 'CellularityPath', 'bcr', 'idc_or_cribriform')
	];
rownames(clinical.filtered) <- substr(rownames(clinical.filtered), 1, 8);

sample.k.order <- names(sample.cluster[[sample.k]]$consensusClass[sample.cluster[[sample.k]]$consensusTree$order])
sample.order <- data.frame(
	samples = sample.k.order,
	clusters = unname(sample.cluster[[sample.k]]$consensusClass[sample.k.order]),
	stringsAsFactors = FALSE
	)
sample.ordered <- sample.order[order(sample.order$clusters), ];
write.table(
    sample.ordered,
    generate.filename('m6A', paste0(input, '_', level, '_sample_clusters'), 'tsv'),
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
    );

gene.k.order <- names(gene.cluster[[gene.k]]$consensusClass[gene.cluster[[gene.k]]$consensusTree$order]);
gene.order <- data.frame(
	genes = gene.k.order,
	clusters <- unname(gene.cluster[[gene.k]]$consensusClass[gene.k.order]),
	stringsAsFactors = FALSE
	);
gene.ordered <- gene.order[order(gene.order$clusters), ];

# this reorders the samples and genes in your matrix based on the consensus clustering results
data <- scale.meth[match(gene.ordered$genes, rownames(scale.meth)), match(sample.ordered$samples, colnames(scale.meth))];
data <- as.data.frame(data);

# create covariates for the heatmap
gene.col <- default.colours(gene.k)[gene.ordered$clusters];
sample.col <- default.colours(sample.k)[sample.ordered$clusters];

sample.covariate <- list(
	rect = list(
		col = 'transparent',
		fill = sample.col,
		lwd = 1.5
		)
	);

gene.covariate <- list(
	rect = list(
		col = 'transparent',
		fill = rev(gene.col),
		lwd = 1.5
		)
	);

# split samples and genes into separate matrices and save in list
gene.split <- split(data, as.factor(gene.ordered$cluster));
sample.and.gene <- lapply(gene.split, function(x) split(as.data.frame(t(x)), sample.ordered$clusters));

# create m6A and covariate heatmaps (create.subset.heatmap and create.covariate.heatmap functions are in file sourced at beginning of script)
plot.objects <- list()
counter <- 1;
for (i in 1:sample.k) {
	for (j in 1:gene.k) {
		plot.objects[[counter]] <- create.subset.heatmap(sample.and.gene[[j]][[i]]);
		counter <- counter + 1;
		}
	annot <- clinical.filtered[match(rownames(sample.and.gene[[1]][[i]]), rownames(clinical.filtered)), ];
	plot.objects[[counter]] <- create.covariate.heatmap(annot);
	counter <- counter + 1;
	}

cov.length <- 500;
gene.length <- nrow(scale.meth) + cov.length;
sample.length <- ncol(scale.meth);

# need to add white space after final sample.k etc is picked
ylab.label <- rev(paste0('P', 1:sample.k));
xlab.label <- paste0('M', 1:gene.k);

# want each row / col to be the same height so scale based on how many samples or gene are in each cluster
panel.heights <- rep(NA, sample.k);
for (i in 1:sample.k) {
	panel.heights[i] <- nrow(sample.and.gene[[1]][[i]]) / sample.length
	}

panel.widths <- rep(NA, gene.k);
for (i in 1:gene.k) {
	panel.widths[i] <- ncol(sample.and.gene[[i]][[1]]) / gene.length
	}

# rev panel heights because boxes are plotted from the bottom
create.multiplot(
	plot.objects = plot.objects,
	file = 'Figure1A.png',
	plot.layout = c(gene.k + 1, sample.k),
	panel.heights = rev(panel.heights),
	panel.widths = c(panel.widths, cov.length / gene.length),
	y.spacing = -0.6,
	x.spacing = 0.2,
	x.relation = 'free',
	y.relation = 'free',
	ylab.label = c('\t    P5', '\t\tP4', '\t    P3', '\tP2', 'P1'),
	xlab.label = c('M1', 'M2', 'M3', 'M4\t    ', 'M5\t\t'),
	xlab.padding = 0,
	xlab.to.xaxis.padding = -1.5,
	ylab.padding = 1,
	xlab.cex = 1,
	ylab.cex = 1,
	yaxis.tck = 0,
	xaxis.tck = 0,
	print.new.legend = TRUE,
	right.padding = 14,
	use.legacy.settings = TRUE,
	legend = list(
		inside = list(
			x = 1.05,
			y = 0.99,
			fun = legend.grob(
				covariates.legend,
				between.row = 0.5,
				size = 1.5,
				label.cex = 0.70,
				title.cex = 0.70,
				title.just = 'left'
				)
			)
		),
	resolution = 200
	);
