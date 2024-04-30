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

clinical.filtered$tcategory <- factor(
	x = substr(
		x = clinical.filtered$pathologic_t,
		start = 1,
		stop = 2
		)
	);
clinical.filtered$isup.grouped <- c(1, 2, 3, 4, 4)[match(clinical.filtered$pathologic_isup_grade, c(1, 2, 3, 4, 5))];

signif <- read.delim(
	file = '2023-08-01_m6A_peaks_with_significant_clinical_associations.txt',
	as.is = TRUE
	);

data.overlap <- data[match(signif$peak.id, rownames(data)), match(patient.overlap, colnames(data))];

id <- 'ENSG00000165671.20'
gene <- 'NSD1'
idc.status <- na.omit(
	ifelse(
		test = clinical.filtered$idc_or_cribriform,
		yes = 'Yes',
		no = 'No'
		)
	);

plot.my.boxplot(
	var1 = idc.status,
	var2 = data.overlap[id, -which(is.na(clinical.filtered$idc_or_cribriform))],
	labels = c('IDC/CA ', paste0(id, '   ', gene)),
	ylab.label = as.expression(bquote(bold(.(gene)~'m'^'6'*'A abundance'))),
	test = 'utest',
	status.one = 'Yes',
	status.two = 'No',
	file = 'Figure2J.png',
	main.cex = 0
	);

id  <- 'ENSG00000196544.8'
gene <- 'BORCS6'

col <- force.colour.scheme(clinical.filtered$pathologic_isup_grade, scheme = 'isup.grade');
col[which(col == 'cornsilk')] <- 'black'

plot.my.boxplot(
	var1 = clinical.filtered$isup.grouped,
	var2 = data.overlap[id, ],
	labels = c('pISUP', paste0(id, '   ', gene)),
	ylab.label = as.expression(bquote(bold(.(gene)~'m'^'6'*'A abundance'))),
	test = 'aov',
	file = 'Figure2I.png',
	col = col,
	main.cex = 0
	);
