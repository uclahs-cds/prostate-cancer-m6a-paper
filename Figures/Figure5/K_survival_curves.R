library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.utilities);

mask <- 'hard';
level <- 'peak';
time.thresh <- '5_year';

file <- 'MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakCounts_HardMask.Adjusted.tsv'

data <- read.delim(
	file = file,
	as.is = TRUE
	);

clinical <- read.delim(
	file = '2020-02-25_CPC-GENE_clinical_from_database.txt',
	as.is = TRUE
	);

clinical.filtered <- clinical[match(colnames(data), clinical$patient_id), c('sample_id', 'bcr', 'time_to_bcr')];

# remove sample without survival data
clinical.survival <- clinical.filtered[-which(is.na(clinical.filtered$time_to_bcr) | is.na(clinical.filtered$bcr)), ];
samples <- substr(clinical.survival$sample_id, 1, 8);

data.survival <- data[, match(samples, colnames(data))];

# remove peaks that aren't seen in more than 5 samples
binary <- ifelse(data.survival == 0, 0, 1);
count <- rowSums(binary);
data.filtered <- data.survival[which(count > 5), ];

df <- ifelse(
	test = data.filtered == 0,
	yes = 0,
	no = 1
	);

surv <- Surv(
	time = clinical.survival$time_to_bcr,
	event = clinical.survival$bcr
	);

cox <- data.frame(
	hr = numeric(),
	lower.95 = numeric(),
	upper.95 = numeric(),
	num.peaks = numeric(),
	pvalue = numeric(),
	stringsAsFactors = FALSE
	);

for (i in 1:nrow(df)) {
	groups <- factor(x = df[i, ], levels = c(0, 1), labels = c(0, 1));

	fit <- fit.coxmodel(
		groups = groups,
		survobj = surv
		);

	cox[i, 'hr'] <- fit[1];
	cox[i, 'lower.95'] <- fit[2];
	cox[i, 'upper.95'] <- fit[3];
	cox[i, 'num.peaks'] <- length(which(groups == 1));
	cox[i, 'pvalue'] <- fit[4];
	}

cox$qvalue <- p.adjust(
	p = cox$pvalue,
	method = 'fdr'
	);

genes <- read.delim(
	file = 'gencode.v34.gene_names.tsv',
	as.is = TRUE,
	header = FALSE
	);

peak.genes <- unlist(
	lapply(
		X = rownames(data.filtered),
		FUN = function(x) {
			strsplit(
				x = x,
				split = ':'
				)[[1]][1]
			}
		)
	);
genes.all <- genes[match(peak.genes, genes$V1), ];
colnames(genes.all) <- c('ensg_id', 'gene.name');

rownames(cox) <- rownames(data.filtered);
cox$gene.name <- genes.all$gene.name;

cox.ordered <- cox[order(cox$pvalue), ];

plot.my.curve <- function(group, surv.obj, peak, gene, level, add.label, file, ...) {

	if (add.label) {
		ylab.label <- 'BCR-Free Rate';
	} else {
		ylab.label <- '';
		}

	obj <- create.km.plot(
		survival.object = surv.obj,
		patient.groups = group,
		resolution = 300,
		statistical.method = 'cox',
		cox.zph.threshold = 0.05,
		ph.assumption.check = 'warning.and.plot',
		show.risktable = FALSE,
		xlab.label = 'Time (Years)',
		ylab.label = ylab.label,
		main.cex = 0.7,
		key.stats.y.pos = 0.15,
		use.legacy.settings = FALSE,
		return.statistics = FALSE,
		risktable.fontsize = 9,
		risk.labels = c('No\nPeak', 'Peak'),
		key.groups.labels = c(paste('No', gene, 'peak'), paste(gene, 'peak')),
		key.groups.title.cex = 0.5,
		key.stats.cex = 0.35,
		yaxis.cex = 0.5,
		xaxis.cex = 0.5,
		xlab.cex = 0.5,
		ylab.cex = 0.5,
		key.groups.cex = 0.4,
		lwd = 0.9,
		censoring.pch.cex = 0.6,
		...
		);

	obj$x.scales$tck <- c(.5, 0);
	obj$y.scales$tck <- c(.5, 0);
	obj$par.settings$axis.components <- list(top = list(pad1 = -2), bottom = list(pad2 = -0.65), left = list(pad2 = -0.65));
	obj$legend$inside$args[[1]]$lines$size = 0.8;
	pdf(file, height = 1.9, width = 2);
	print(obj);
	dev.off();
	}

# Figure 5K
VCAN.peak <- 'ENSG00000038427.16:9';
plot.my.curve(
	group = as.numeric(df[match(VCAN.peak, rownames(df)),]),
	file = generate.filename('m6A', paste0('VCAN_kaplan_meier_curve_5_year_', time.thresh == '5_year'), 'pdf'),
	surv.obj = surv,
	peak = VCAN.peak,
	gene = 'VCAN',
	level = 'peak',
	add.label = TRUE,
	key.stats.x.pos = .61,
	xlimits = if (time.thresh == '5_year') c(0, 5) else NA,
	xat = if (time.thresh == '5_year') seq(0, 5, 1) else NA,
	predefined.p = if (time.thresh == '5_year') cox[VCAN.peak, ]$pvalue else NULL,
	predefined.hr = if (time.thresh == '5_year') cox[VCAN.peak, ]$hr else NULL,
	predefined.hr.ci = if (time.thresh == '5_year') as.numeric(cox[VCAN.peak, c('lower.95', 'upper.95')]) else NULL
	);

# Supplementary Figure 5G
INHBA.peak <- 'ENSG00000122641.11:1';
plot.my.curve(
	group = as.numeric(df[match(INHBA.peak, rownames(df)),]),
	file = generate.filename('m6A', paste0('INHBA_kaplan_meier_curve_5_year_', time.thresh == '5_year'), 'pdf'),
	surv.obj = surv,
	peak = INHBA.peak,
	gene = 'INHBA',
	level = 'peak',
	add.label = TRUE,
	key.stats.x.pos = .6,
	xlimits = if (time.thresh == '5_year') c(0, 5) else NA,
	xat = if (time.thresh == '5_year') seq(0, 5, 1) else NA,
	predefined.p = if (time.thresh == '5_year') cox[INHBA.peak, ]$pvalue else NULL,
	predefined.hr = if (time.thresh == '5_year') cox[INHBA.peak, ]$hr else NULL,
	predefined.hr.ci = if (time.thresh == '5_year') as.numeric(cox[INHBA.peak, c('lower.95', 'upper.95')]) else NULL
	);

# Supplementary Figure 5H
ZFHX4.peak <- 'ENSG00000091656.19:8';
plot.my.curve(
	group = as.numeric(df[match(ZFHX4.peak, rownames(df)),]),
	file = generate.filename('m6A', paste0('ZFHX4_kaplan_meier_curve_5_year_', time.thresh == '5_year'), 'pdf'),
	surv.obj = surv,
	peak = ZFHX4.peak,
	gene = 'ZFHX4',
	level = 'peak',
	add.label = TRUE,
	key.stats.x.pos = .56,
	xlimits = if (time.thresh == '5_year') c(0, 5) else NA,
	xat = if (time.thresh == '5_year') seq(0, 5, 1) else NA,
	predefined.p = if (time.thresh == '5_year') cox[ZFHX4.peak, ]$pvalue else NULL,
	predefined.hr = if (time.thresh == '5_year') cox[ZFHX4.peak, ]$hr else NULL,
	predefined.hr.ci = if (time.thresh == '5_year') as.numeric(cox[ZFHX4.peak, c('lower.95', 'upper.95')]) else NULL
	);
