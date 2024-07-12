library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

source('plot.general.contingency.table.R');

level <- 'peak'
input <- 'adjusted'
sample.k <- 5
adjust <- FALSE

sample.part <- paste0(input, '_', level, '_sample_ConsensusClusterPlus');
sample.file <- paste0(sample.part, '/', sample.part, '.k=', sample.k, '.consensusClass.csv');
sample.clusters <- read.delim(
	file = sample.file,
	as.is = TRUE,
	header = FALSE,
	sep = ','
	);

clinical <- read.delim(
	file = '2020-02-25_CPC-GENE_clinical_from_database.txt',
	as.is = TRUE
	);

annot <- clinical[match(sample.clusters$V1, clinical$patient_id), ];
annot$cluster <- paste0('P', sample.clusters$V2);
annot$isup.grouped <- c(1, 2, 3, 4, 4)[match(annot$pathologic_gleason_grades, c('3+3', '3+4', '4+3', '4+4', '4+5'))];

covs.to.keep <- c('pISUP', 'PSA', 'pT', 'Age', 'IDC/CA');
pvalue <- rep(NA, length(covs.to.keep));

# ISUP - two categorical variables
pvalue[1] <- chisq.test(
	x = table(
		sample.clusters$V2,
		annot$isup.grouped
		),
	simulate.p.value = FALSE
	)$p.value;

# PSA - continuous vs categorical
pvalue[2] <- summary(
	aov(
		formula = sample.clusters$V2 ~ annot$pre_treatment_psa
		)
	)[[1]][['Pr(>F)']][1];

# Tcat - two categorical variables
pvalue[3] <- chisq.test(
	x = table(
		sample.clusters$V2,
		substr(annot$pathologic_t, 1, 2)
		)
	)$p.value;

# Age - continuous vs categorical
pvalue[4] <- summary(
	aov(
		formula = sample.clusters$V2 ~ annot$age_at_treatment
		)
	)[[1]][['Pr(>F)']][1];

# IDC - two categorical variables
pvalue[5] <- chisq.test(
	x = table(
		sample.clusters$V2,
		annot$idc_or_cribriform
		)
	)$p.value

pvals.adj <- pvalue;
xlab.label <- 'P value\n';
ord <- order(pvals.adj);

pvals.clin <- data.frame(
	pvalue = -log10(pvals.adj[ord]),
	order = 1:length(covs.to.keep),
	stringsAsFactors = FALSE
	);

cluster.barplot <- create.barplot(
	formula = rev(pvalue) ~ order,
	data = pvals.clin,
	resolution = 500,
	xaxis.tck = c(0.2, 0),
	yaxis.tck = c(0.2, 0),
	ylab.cex = 0.35,
	ylab.label = xlab.label,
	ylab.axis.padding = -1.5,
	yaxis.lab = c(expression(bold('10')^bold('0')), expression(bold('10')^bold('-1')), expression(bold('10')^bold('-2'))),
	xlab.cex = 0,
	xaxis.lab = rev(covs.to.keep[ord]),
	xaxis.rot = 0,
	yaxis.cex = 0.275,
	xaxis.cex = 0.275,
	yat = seq(0, 2, 1),
	ylimits = c(0, 2),
	abline.h = 1.30103,
	abline.col = 'red',
	bottom.padding = 0,
	top.padding = 0,
	plot.horizontal = FALSE,
	axes.lwd = .25
	);
cluster.barplot$par.settings$axis.components <- list(
	bottom = list(pad1 = 0.3),
	left = list(pad1 = 0.3)
	);
png(
	filename = Figure2B.png',
	height = 1,
	width = 2.25,
	units = 'in',
	res = 500
	)
print(cluster.barplot)
dev.off()

plot.data <- table(
	paste0('P', sample.clusters$V2),
	ifelse(
		test = annot$idc_or_cribriform,
		yes = 'Yes',
		no = 'No'
		)
	);

IDC.plot <- plot.general.contingency.table(
	plot.data = plot.data,
	filename = 'Figure2D.png',
	height = 2,
	width = 1.65,
	xlabel = 'IDC/CA',
	ylabel = '',
	left.padding = 1.5,
	y.spacing = -0.2,
	x.spacing = -0.2,
	panel.widths = c(1, 0.3),
	main.x = 0.4,
	main.cex = 0.45,
	text.cex = 0.5,
	xlab.cex = 0.4,
	ylab.cex = 0.45,
	xaxis.cex = 0.4,
	yaxis.cex = 0.4
	);

plot.data <- table(
	paste0('P', sample.clusters$V2),
	substr(annot$pathologic_t, 1, 2)
	);

plot.general.contingency.table(
	plot.data = plot.data,
	filename = 'Figure2C.png',
	xlabel = 'pT',
	ylabel = '',
	height = 2,
	width = 1.65,
	left.padding = 1.5,
	y.spacing = -0.2,
	x.spacing = -0.1,
	panel.widths = c(1, 0.3),
	main.x = 0.4,
	main.cex = 0.45,
	text.cex = 0.5,
	xlab.cex = 0.4,
	ylab.cex = 0.45,
	xaxis.cex = 0.4,
	yaxis.cex = 0.4
	);
