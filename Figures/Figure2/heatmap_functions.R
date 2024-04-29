create.subset.heatmap <- function(data) {

	key.min <- -3;
	key.max <- 3;
	key.colour.interval.num <- 50;
	key.scale <- c(
		seq(key.min, 0, -key.min / key.colour.interval.num),
		seq(0, key.max, key.max / key.colour.interval.num)
		);
	key.scale <- unique(key.scale);

	rnaseq.heatmap <- create.heatmap(
		x = data,
		clustering.method = 'none',
		at = key.scale,
		colour.scheme = c('darkorchid4', 'white', 'darkgreen'),
		print.colour.key = FALSE,
		same.as.matrix = TRUE,
		scale.data = FALSE
		);

	return(rnaseq.heatmap)
	}

create.covariate.heatmap <- function(annot) {
	numeric <- annot;

	numeric['age_at_treatment'] <- force.colour.scheme(annot$age_at_treatment, scheme = 'age.categorical.prostate');
	numeric['pre_treatment_psa'] <- force.colour.scheme(annot$pre_treatment_psa, scheme = 'psa.categorical');
	numeric['clinicalt3'] <- force.colour.scheme(annot$clinicalt3, scheme = 'clinicalt3');
	numeric['pathologic_isup_grade'] <- force.colour.scheme(annot$pathologic_isup_grade, scheme = 'isup.grade');

	numeric[which(annot$bcr == TRUE), 'bcr'] <- 'black'
	numeric[which(annot$bcr == FALSE), 'bcr'] <- 'white'

	numeric[which(annot$idc_or_cribriform == TRUE), 'idc_or_cribriform'] <- 'darkred'
	numeric[which(annot$idc_or_cribriform == FALSE), 'idc_or_cribriform'] <- 'white'

        numeric[which(annot$CellularityPath == 0.5), 'CellularityPath'] <- '#bcf5f9'
        numeric[which(annot$CellularityPath == 0.6), 'CellularityPath'] <- '#89c5fd'
        numeric[which(annot$CellularityPath == 0.7), 'CellularityPath'] <- '#3a80ec'
        numeric[which(annot$CellularityPath == 0.8), 'CellularityPath'] <- '#0229bf'
        numeric[which(annot$CellularityPath >= 0.9), 'CellularityPath'] <- '#080b6c'

	heatmap.colours <- na.omit(unique(unlist(numeric)));
	for (i in 1:length(heatmap.colours)) {
		numeric[numeric == heatmap.colours[i]] <- i;
		}
	numeric <- sapply(numeric, as.numeric);

	cov.heatmap <- create.heatmap(
		x = numeric,
		clustering.method = 'none',
		same.as.matrix = TRUE,
		colour.scheme = heatmap.colours,
		grid.col = TRUE,
		grid.row = FALSE,
		total.colours = length(heatmap.colours) + 1,
		fill.colour = 'darkgray',
		print.colour.key = FALSE
		);

	return(cov.heatmap);
	}

# set general variables
key.min <- -3;
key.max <-  3;
key.colour.interval.num <- 50;
key.scale <- c(
	seq(key.min, 0, -key.min / key.colour.interval.num),
	seq(0, key.max, key.max / key.colour.interval.num)
	);
key.scale <- unique(key.scale);

covariates.legend <- list(
	legend = list(
		colours = c('gray75', 'gray50', 'gray25', 'gray0'),
		labels = c('40 - 50', '50 - 65', '65 - 70', expression(''>= 70)), #nolint
		title = expression(bold(underline('Age')))
		),
	legend = list(
		colours = c('#FEE6CE', '#FDAE6B'),
		labels = c('0 - 9.9', '10 - 19.9'),
		title = expression(bold(underline('PSA (ng/mL)')))
		),
	legend = list(
		colours = c('#6DC46E', '#2F6D60'),
		labels = c('T2', 'T3'),
		title = expression(bold(underline('pT')))
		),
	legend = list(
		colours = c('cornsilk', 'yellow', 'orange', 'maroon3', 'red'),
		labels = c('1', '2', '3', '4', '5'),
		title = expression(bold(underline('pISUP')))
		),
	legend = list(
		colours = c('#bcf5f9', '#89c5fd', '#3a80ec', '#0229bf', '#080b6c'),
		labels = c('0.5', '0.6', '0.7', '0.8', expression(''>= 0.9)),
		title = expression(bold(underline('Cellularity')))
		),
	legend = list(
		colours = c('black', 'white'),
		labels = c('Yes', 'No'),
		title = expression(bold(underline('BCR')))
		),
	legend = list(
		colours = c('darkred', 'white'),
		labels = c('Yes', 'No'),
		title = expression(bold(underline('IDC/CA')))
		),
	legend = list(
		colours = c('darkgreen', 'white', 'darkorchid4'),
		labels = c(
			as.expression(bquote(''<=.(round(key.scale[1], digits = 2)))),
			as.character(median(key.scale)),
			as.expression(bquote(''>=.(round(key.scale[length(key.scale)], digits = 2))))
			),
		title = expression(bold(underline('Z'[Intensity]))),
		continuous = TRUE,
		tck = 0,
		at = c(0, 50, 100)
		)
	);
