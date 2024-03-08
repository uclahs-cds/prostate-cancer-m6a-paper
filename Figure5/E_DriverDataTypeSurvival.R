library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

config.path <- '/hot/users/rhughwhite/projects/PRAD-000080-N6Methyl/code/project-ProstateCancer-m6A/survival/config/drivers_data_type_comparison_forest_plot_continuous_m6A.config';
config <- read.config.file(config.path);
survival.results <- read.delim(config$survival.results);
max.gene.results <- aggregate(survival.results$hr, by = list(gene = survival.results$gene), FUN = max, na.rm = TRUE);
survival.results$max.gene.hr <- max.gene.results$x[match(survival.results$gene, max.gene.results$gene)];
survival.results$data.type <- sub('.binary', '', survival.results$data.type);
survival.results$data.type <- factor(survival.results$data.type, levels = c('Gain', 'Loss', 'RNA', 'm6A', 'protein'));
survival.results <- survival.results[order(survival.results$data.type, survival.results$hr, decreasing = FALSE), ];

output <- paste0(config$output.path, '_binary_m6A_', as.logical(config$binary.m6A), '_drivers_data_type_comparison_forest_plot.pdf');

survival.results$label <- survival.results$gene;
survival.results$log.HR <- log2(survival.results$hr);
survival.results$log.HR.high <- log2(survival.results$upper.95);
survival.results$log.HR.low <- log2(survival.results$lower.95);
survival.results$order <- 1:nrow(survival.results);
survival.results <- survival.results[!is.na(survival.results$hr), ];

x.at.values <- seq(-4, 4, 2);
x.at.labels <- sapply(x.at.values, function(power) as.expression(bquote(bold(.('2') ^.(power)))));
forest <- create.scatterplot(
	formula = 1:nrow(survival.results) ~ log.HR,
	data = survival.results,
	ylab.label = 'Driver',
	xlab.label = 'Hazard Ratio',
    yaxis.lab = NULL,
    ylimits = c(0.5, nrow(survival.results) + .5), #need to add .5 spacing to line up with barplot
	xlimits = log2(c(0.05, 2^5)),
    xaxis.lab = x.at.labels,
    xat = log2(sapply(x.at.values, function(x) 2^x)),
	xlab.cex = .65,
    xaxis.cex = .6,
	yaxis.cex = .4,
	yaxis.tck = 0,
    xaxis.tck = c(.5, 0),
	ylab.cex = 0,
	abline.v = 0,
	abline.col = 'grey',
	abline.lwd = 2,
	abline.lty = 3,
	pch = 19,
	#fill = 'transparent',
	error.bar.lwd = 1,
	error.whisker.angle = 0,
	x.error.right = abs(survival.results$log.HR - survival.results$log.HR.high),
	x.error.left = abs(survival.results$log.HR - survival.results$log.HR.low),
	y.error.bar.col = 'black'
	);

data.type.colours <- read.delim('input_data/data_types_colours.txt');
data.type.colours <- setNames(data.type.colours$data.type.colours, row.names(data.type.colours));
data.type.colours['Gain'] <- 'firebrick1';
data.type.colours['Loss'] <- 'dodgerblue';
data.type.colours <- data.type.colours[as.character(levels(survival.results$data.type))];

data.type.legend.labels <- names(data.type.colours);
data.type.legend.labels[data.type.legend.labels == 'm6A'] <- expression(m^6 * A);

data.types.map <- create.heatmap(
    x = data.frame(colour = data.type.colours[as.character(survival.results$data.type)]),
    clustering.method = 'none',
    yaxis.lab = survival.results$label,
    yat = 1:nrow(survival.results),
    yaxis.cex = .5,
    yaxis.tck = 0,
    xaxis.tck = 0,
    same.as.matrix = FALSE,
    print.colour.key = FALSE,
    input.colours = TRUE,
    axes.lwd = 1,
    );

data.types.legend <- list(
    colours = rev(data.type.colours),
    labels = rev(data.type.legend.labels),
    title = bquote(bold(underline('Data type'))),
    cex = .7
    );

legend <- legend.grob(
    legends = list(legend = data.types.legend),
    label.cex = 0.5,
    title.cex = 0.6,
    title.just = 'left',
    title.fontface = 'plain',
    size = 1.5
    );

FDR.barplot <- FDR.barplot(
    FDR = survival.results$pvalue,
    FDR.threshold = 0.05,
    log10 = TRUE,
    FDR.intervals = -log10(c(.1, 0.001, 0.00001)),
    plot.xaxis.labels = TRUE,
    xlab.label = 'P',
    xlab.cex = .65,
    xaxis.cex = .6
    );

plot.objects.widths <- c(0.2, 0.6, 0.30);

create.multipanelplot(
	plot.objects = list(data.types.map, forest, FDR.barplot),
	file = output,
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = plot.objects.widths,
	x.spacing = c(-.8, -.5),
	resolution = 400,
	legend = list(
		right = list(
			fun = legend
			)
		),
	left.padding = -2,
    right.padding = 0,
	right.legend.padding = 0,
    top.legend.padding = 0,
	height = 4,
	width = 4
	);
