library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(scales);

rna <- read.delim(
	file = '2021-11-16_RSEM_gene_TPM.txt',
	as.is = TRUE,
	check.names = FALSE
	);

protein <- read.delim(
	file = '2018-04-06_prot_combat_imputed.txt',
	as.is = TRUE
	);

peaks <- read.delim(
	file = 'MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakSum_GeneCounts_Soft.Adjusted.tsv',
	as.is = TRUE
	);
peaks.binary <- ifelse(peaks > 0, 1, 0);
peak.count.per.gene <- data.frame(
	peak.count = rowSums(peaks.binary),
	stringsAsFactors = FALSE
	);

methyl <- read.delim(
	file = 'MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakSum_GeneCounts.Adjusted.tsv',
	as.is = TRUE
	);

methyl.per.gene <- data.frame(
	mean = apply(methyl, 1, mean),
	median = apply(methyl, 1, median),
	stringsAsFactors = FALSE
	);

# protein & RNA correlation
gene.overlap <- intersect(protein$gene.names, rna$Symbol);
sample.overlap <- intersect(colnames(protein), substr(colnames(rna), 1,8));

protein.filtered <- protein[match(gene.overlap, protein$gene.names), match(sample.overlap, colnames(protein))];
rna.filtered <- rna[match(gene.overlap, rna$Symbol), match(sample.overlap, substr(colnames(rna), 1, 8))];

results <- data.frame(
	gene = gene.overlap,
	rho = rep(NA, length = nrow(protein.filtered))
	);

for (i in 1:nrow(protein.filtered)) {
	results[i, 'rho'] <- cor(
		x = as.numeric(protein.filtered[i,]),
		y = as.numeric(rna.filtered[i,]),
		method = 'spearman'
		);
	}

gene.annot <- read.delim(
	file = 'gencode.v34.gene_names.tsv',
	as.is = TRUE,
	header = FALSE
	);
	
peak.count.per.gene$gene <- gene.annot$V2[match(rownames(peak.count.per.gene), gene.annot$V1)];
methyl.per.gene$gene <- gene.annot$V2[match(rownames(methyl.per.gene), gene.annot$V1)];
peak.overlap <- intersect(results$gene, peak.count.per.gene$gene);

results.filtered <- results[match(peak.overlap, results$gene), ];
peak.count.per.gene.filtered <- peak.count.per.gene[match(peak.overlap, peak.count.per.gene$gene), ];
methyl.per.gene.filtered <- methyl.per.gene[match(peak.overlap, methyl.per.gene$gene), ];

toplot <- data.frame(
	x = peak.count.per.gene.filtered$peak.count,
	y = results.filtered$rho
	);

create.hexbinplot(
	formula = y ~ x,
	data = toplot,
	mincnt = 1,
	filename = 'Supplementary_Figure2M.png',
	ylab.label = bquote(bold("Spearman's"~rho~"(RNA-Protein)")),
	xlab.label = expression(bold('# of samples with a peak')),
	xaxis.cex = 0.7,
	yaxis.cex = 0.7,
	ylab.cex = 1,
	xlab.cex = 1,
	xaxis.tck = 0.2,
	yaxis.tck = 0.2,
	ylimit = c(-1, 1),
	yat = seq(-1, 1, 0.5),
	xlimit = c(-2, 150),
	xat = seq(0, 150, 50),
	aspect = 1,
	key = get.corr.key(
		x = toplot$x,
		y = toplot$y,
		label.items = c('spearman', 'spearman.p'),
		key.cex = 1,
		x.pos = 0.01,
		y.pos = 0.05
		),
	right.padding = 1.5,
	width = 5,
	height = 5,
	size.units = 'in',
	resolution = 300
	);

num.of.splits <- 10
n.peaks <- length(peak.count.per.gene.filtered$peak.count);
decile <- rep(
        x = 1:num.of.splits,
        each = (n.peaks / num.of.splits)
        );
decile[length(decile):n.peaks] <- num.of.splits;
peak.count.ordered <- peak.count.per.gene.filtered[order(peak.count.per.gene.filtered$peak.count, decreasing = TRUE), ];
rho.ordered <- results.filtered[match(peak.count.ordered$gene, results.filtered$gene), ];

rank <- data.frame(
        rho = rho.ordered$rho,
        decile = decile
        );

create.boxplot(
        formula = rho ~ as.factor(decile),
        data = rank,
        file = 'Supplmentary_Figure2K.png',
        add.stripplot = TRUE,
        xaxis.cex = 1,
        yaxis.cex = 1,
	yat = seq(-0.5, 1, 0.5),
	ylimit = c(-0.5, 1),
        ylab.label = bquote(bold("Spearman's"~rho~"(RNA-Protein)")),
        ylab.cex = 1,
        xlab.label = 'Decile\n(ranked by # of samples with a peak)',
        xlab.cex = 1,
	resolution = 200
        );

toplot.mean <- data.frame(
	x = log2(methyl.per.gene.filtered$mean),
	y = results.filtered$rho,
	z = log2(methyl.per.gene.filtered$median)
	);

create.hexbinplot(
	formula = y ~ x,
	data = toplot.mean,
	maxcnt = 50,
	mincnt = 1,
	filename = 'Supplementary_Figure2L.png',
	ylab.label = bquote(bold("Spearman's"~rho~"(RNA-Protein)")),
	xlab.label = expression(bold('Mean log')[bold('2')] * bold('(m')^bold('6') * bold('A abundance)')),
	xaxis.fontface = 'bold',
	yaxis.fontface = 'bold',
	xaxis.cex = 0.7,
	yaxis.cex = 0.7,
	ylab.cex = 1,
	xlab.cex = 1,
	xaxis.tck = 0.2,
	yaxis.tck = 0.2,
	ylimit = c(-1, 1),
	yat = seq(-1, 1, 0.5),
	xlimit = c(-5, 20),
	xat = seq(-5, 20, 5),
	aspect = 1,
	key = get.corr.key(
		x = toplot.mean$x,
		y = toplot.mean$y,
		label.items = c('spearman', 'spearman.p'),
		key.cex = 1,
		x.pos = 0.01,
		y.pos = 0.05
		),
	right.padding = 1.5,
	width = 5,
	height = 5,
	size.units = 'in',
	resolution = 300
	);
