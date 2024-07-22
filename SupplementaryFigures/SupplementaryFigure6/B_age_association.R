library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

source('zoom_in_plot_functions.R');

data <- read.delim(
	file = 'MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakCounts.Adjusted.tsv',
	as.is = TRUE
	);

clinical <- read.delim(
	file = '2020-02-25_CPC-GENE_clinical_from_database.txt',
	as.is = TRUE
	);

patient.overlap <- intersect(colnames(data), clinical$patient_id);
cols.to.keep  <- c('pga', 'idc_or_cribriform', 'pathologic_isup_grade', 'pathologic_t', 'age_at_treatment', 'pre_treatment_psa');
clinical.filtered <- clinical[match(patient.overlap, clinical$patient_id), cols.to.keep];

signif <- read.delim(
	file = '2023-08-01_m6A_peaks_with_significant_clinical_associations.txt',
	as.is = TRUE
	);

data.overlap <- data[match(signif$peak.id, rownames(data)), match(patient.overlap, colnames(data))];

gene <- 'TAOK2'
id <- 'ENSG00000149930.18:7'
xat <- seq(0, 300, 100);
xaxis.lab <- xat;
key.y <- .2;
plot.my.scatterplot(
	var1 = clinical.filtered$age_at_treatment,
	var2 = data.overlap[id, ],
	labels = c('Age (at treatment)', paste0(id, '   ', gene)),
	file = 'SupplementaryFigure6B.png',
	xat = xat,
	xaxis.lab = xaxis.lab,
	xlab.label = as.expression(bquote(bold(.(gene)~'m'^'6'*'A abundance'))),
	key.y = key.y,
	main.cex = 0,
	height = 3,
	width = 2.5
	);
