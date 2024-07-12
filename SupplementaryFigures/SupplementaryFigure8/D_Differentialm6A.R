library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(stringr);
library(biomaRt);
# general parameters
script.name <- '03_driver_peak_cluster_heatmap_nonclonal.R';
timestamp <- format(Sys.time(), '%Y%m%d-%H%M%S');

config.path <- 'code/project-ProstateCancer-m6A/somatic_driver_events/configs/03_driver_peak_cluster_heatmap.config';
config <- read.config.file(config.path);
output.dir <- config$output.path;

if (!dir.exists(config$output.path)) {
	dir.create(config$output.path, recursive = TRUE);
	}
session.info.dir <- paste0(output.dir, '/session_info/');

### ANALYSIS ###########################################################

peak.clusters <- read.delim(config$peak.clusters);
row.names(peak.clusters) <- peak.clusters$peak;

driver.clusters <- read.delim(config$driver.clusters);
driver.clusters$driver <- sub('ETS', 'Loss.ETS', driver.clusters$driver);
row.names(driver.clusters) <- driver.clusters$driver;

peak.K <- max(peak.clusters$cluster);
plot.path <- paste0(config$output.path, timestamp, '_driver_peak_K_', peak.K, '_cluster_heatmap.pdf');

beta <- read.delim(config$beta);
FDR <- read.delim(config$FDR);

names(beta) <- sub('ETS', 'Loss.ETS', names(beta));
names(FDR) <- sub('ETS', 'Loss.ETS', names(FDR));

#define colourkey for scaled MI
beta.colours <- c('#2166ac', '#59a8f7', '#b9d9fa', 'white', '#fcc2c9', '#e65365', '#b2182b');
#Use a colour scheme with set intervals in order to highlight smaller fold changes/betas
beta.breakpoints <- c(-1, -0.1, 0, 0.1, 1);
beta.colours.matrix <- matrix('slategrey', nrow = nrow(beta), ncol = ncol(beta));
beta.colours.matrix[beta >= beta.breakpoints[length(beta.breakpoints)]] <- beta.colours[length(beta.colours)];
beta.colours.matrix[beta >= beta.breakpoints[length(beta.breakpoints) - 1] & beta < beta.breakpoints[length(beta.breakpoints)]] <- beta.colours[length(beta.colours) - 1];
beta.colours.matrix[beta >= beta.breakpoints[length(beta.breakpoints) - 2] & beta < beta.breakpoints[length(beta.breakpoints) - 1]] <- beta.colours[length(beta.colours) - 2];

beta.colours.matrix[beta <= beta.breakpoints[1]] <- beta.colours[1];
beta.colours.matrix[beta <= beta.breakpoints[2] & beta > beta.breakpoints[1]] <- beta.colours[2];
beta.colours.matrix[beta <= beta.breakpoints[3] & beta > beta.breakpoints[2]] <- beta.colours[3];
dimnames(beta.colours.matrix) <- dimnames(beta);
#create corresponding legend
beta.legend <- list(
	colours = rev(beta.colours[which(beta.colours != 'white')]),
	labels = rev(c(
		as.expression(bquote(paste('' <= .(beta.breakpoints[1])))),
		as.expression(bquote(paste('' <= .(beta.breakpoints[2])))),
		as.expression(bquote(paste('' <= .(beta.breakpoints[3])))),
		as.expression(bquote(paste('' >= .(beta.breakpoints[length(beta.breakpoints) - 2])))),
		as.expression(bquote(paste('' >= .(beta.breakpoints[length(beta.breakpoints) - 1])))),
		as.expression(bquote(paste('' >= .(beta.breakpoints[length(beta.breakpoints)]))))
		)),
	title = bquote(underline(bold('log'[2] ~ 'Fold Change'))),
	cex = .7
	);

#Function to create the main heatmap showing mutual information gene cluster
#cluster: which cluster to plot. Used to extract data from gene.cluster and for labelling plot
#gene.clusters: data.frame with gene and cluster columns
#pathways: whether to add pathway plot to the 'right' of cluster plot
#gene.cluster.pathways: data.frame containing significantly enriched pathways in 'source' column
#gene.to.pathway: data.frame with genes as row.names, columns as pathways, logical (TRUE,FALSE) or binary (0,1) indicating if gene is tagged with each pathway term
create.cluster.heatmap <- function(
	peak.cluster,
	peaks,
	drivers,
	beta.df,
	beta.colours.df,
	beta.colours,
	beta.breakpoints,
	print.xlabels,
	xaxis.lab = rep('', length(peaks))
	) {

	#extract cluster-specific data
	#perform peak-beta clustering
	beta.cluster <- beta.df[peaks, rev(drivers), drop = FALSE];
	beta.colours.sub <- beta.colours.df[peaks, rev(drivers), drop = FALSE];
	if (length(peaks) == 1) {
		beta.cluster.order <- data.frame(t(beta.cluster));
		} else {
		beta.dist <- dist(beta.cluster);
		beta.dist[is.na(beta.dist)] <- 0; #some peaks may not have dist mesaures due to missing data - set to 0
		beta.cluster.order <- beta.colours.sub[hclust(beta.dist)$order, ];
		}
	if (print.xlabels) {
		cluster.heatmap <- create.heatmap(
			x = t(beta.cluster.order),
			clustering.method = 'none',
			same.as.matrix = TRUE,
			print.colour.key = FALSE,
			input.colours = TRUE,
			xlab.label = bquote(DMPC * .(peak.cluster)), #dysregulated m6A subgroups
			xlab.cex = .6,
			xaxis.lab = xaxis.lab,
			xaxis.cex = .4,
			xaxis.rot = 90,
			xaxis.fontface = 'bold',
			yaxis.tck = 0,
			xaxis.tck = 0,
			axes.lwd = 1
			);
		} else {
		cluster.heatmap <- create.heatmap(
			x = t(beta.cluster.order),
			clustering.method = 'none',
			same.as.matrix = TRUE,
			print.colour.key = FALSE,
			input.colours = TRUE,
			xlab.cex = 0,
			yaxis.tck = 0,
			xaxis.tck = 0,
			axes.lwd = 1
			);
		}
	return(list(cluster.heatmap)); #returning a list allows easier combination of plots later
	}

num.sig.peaks <- apply(FDR, 2, function(driver) length(which(driver < config$FDR.threshold)));
num.sig.peaks.df <- data.frame(
	driver = names(num.sig.peaks),
	driver.sub = sub('.reader|.writer|.eraser', '', names(num.sig.peaks)),
	num.peaks = num.sig.peaks,
	cluster = driver.clusters[names(num.sig.peaks), 'cluster'],
	row.names = 
	);
num.sig.peaks.df <- num.sig.peaks.df[order(num.sig.peaks.df$cluster, num.sig.peaks.df$num.peaks, decreasing = FALSE), ];

#generate gene cluster heatmaps and pathway heatmaps
cluster.heatmaps <- lapply(
	unique(driver.clusters$cluster),
	function(driver.cluster)
		lapply(
			unique(peak.clusters$cluster),
			function(peak.cluster) {
				peaks <- peak.clusters[peak.clusters$cluster == peak.cluster, 'peak'];
				if (as.logical(config$print.pca.genes) | as.logical(config$print.driver.genes)) {
					xaxis.lab <- peak.clusters[peaks, 'gene'];
					} else {
					xaxis.lab <- rep('', length(peaks));
					}
				create.cluster.heatmap(
					peak.cluster = peak.cluster,
					peaks = peaks,
					drivers = num.sig.peaks.df[num.sig.peaks.df$cluster == driver.cluster, 'driver'],
					beta.df = beta,
					beta.colours.df = beta.colours.matrix,
					beta.colours,
					beta.breakpoints,
					print.xlabels = driver.cluster == max(unique(driver.clusters$cluster)),
					xaxis.lab = xaxis.lab
					);
				}
			)
	);

mutation.legend <- create.mutation.heatmap(mutations = driver.clusters$driver, plot.region = FALSE, clonal = FALSE, fontface = 'bold')$legend;

#treating ETS as a loss so no SVs
mutation.legend$legend$colours <- mutation.legend$legend$colours[names(mutation.legend$legend$colours) != 'SV'];
mutation.legend$legend$labels <- mutation.legend$legend$labels[mutation.legend$legend$labels != 'SV'];

mutation.plots <- lapply(
	unique(driver.clusters$cluster),
	function(driver.cluster)
		create.mutation.heatmap(
			mutations = num.sig.peaks.df[num.sig.peaks.df$cluster == driver.cluster, 'driver.sub'],
			plot.region = FALSE,
			add.labels = TRUE,
			return.legend = FALSE,
			transpose.matrix = FALSE,
			label.cex = .6,
			axes.lwd = 1,
			clonal = FALSE,
			fontface = 'bold'
			)
	);

m6A.enzyme.colours <- setNames(c(default.colours(3), 'slategrey'), c('eraser', 'reader', 'writer', 'N/A'));

#get mutation colour and label mappings
m6A.enzyme.legend <- list(
	colours = unlist(m6A.enzyme.colours),
	labels = str_to_title(names(m6A.enzyme.colours)),
	title = expression(underline(bold(m^6 * A~Enzyme))),
	lwd = 0.5
	);

create.m6A.enzyme.heatmap <- function(mutations, m6A.enzyme.colours) {
	mutation.colours <- rep('slategrey', length(mutations));
	for (enzyme.type in names(m6A.enzyme.colours)) {
		mutation.colours[grep(enzyme.type, mutations)] <- m6A.enzyme.colours[[enzyme.type]];
		}
	mutation.colours <- matrix(mutation.colours);
	mutation.heatmap <- create.heatmap(
		mutation.colours,
		input.colours = TRUE,
		clustering.method = 'none',
		print.colour.key = FALSE,
		xaxis.cex = 0,
		ylab.cex = 0,
		xaxis.tck = 0,
		yaxis.tck = 0,
		axes.lwd = 1
		);
	return(mutation.heatmap);
	}

m6A.enzyme.plots <- lapply(
	unique(driver.clusters$cluster),
	function(driver.cluster)
		create.m6A.enzyme.heatmap(
			mutations = num.sig.peaks.df[num.sig.peaks.df$cluster == driver.cluster, 'driver'],
			m6A.enzyme.colours = m6A.enzyme.colours
			)
	);

xlabels.num.peaks <- c(0, 10^seq(1:floor(log10(max(num.sig.peaks.df$num.peaks)))));
x.at.num.peaks <- log10(xlabels.num.peaks + 1);
num.sig.peaks.df$num.peaks.log <- log10(num.sig.peaks.df$num.peaks + 1);

create.n.peaks.plot <- function(
	drivers,
	num.sig.peaks.df,
	xlabels.num.peaks,
	x.at.num.peaks,
	print.xlabels,
	xlimits
	) {
	num.sig.peaks.sub <- num.sig.peaks.df[num.sig.peaks.df$driver %in% drivers, ];
	if (print.xlabels) {
		driver.n.peaks.plot <- create.barplot(
			1:nrow(num.sig.peaks.sub) ~ num.peaks.log,
			data = num.sig.peaks.sub,
			plot.horizontal = TRUE,
			xlimits = xlimits,
			border.lwd = 0,
			xaxis.tck = 0,
			xaxis.cex = 0.5,
			yaxis.cex = 0,
			yaxis.tck = 0,
			xaxis.lab = xlabels.num.peaks,
			xat = x.at.num.peaks,
			xaxis.rot = 90,
			ylab.label = '',
			left.padding = 0,
			xlab.label = expression(bold('\n\tNumber of peaks\n\t   (FDR < 0.1)')),
			xlab.cex = .6,
			xaxis.fontface = 'bold',
			ylab.cex = 0,
			);
		} else {
		driver.n.peaks.plot <- create.barplot(
			1:nrow(num.sig.peaks.sub) ~ num.peaks.log,
			data = num.sig.peaks.sub,
			plot.horizontal = TRUE,
			xlimits = xlimits,
			border.lwd = 0,
			xaxis.tck = 0,
			xaxis.cex = 0,
			yaxis.cex = 0,
			yaxis.tck = 0,
			ylab.label = '',
			xlab.axis.padding = 0,
			left.padding = 0,
			xlab.cex = 0,
			ylab.cex = 0,
			);
		}
	return(driver.n.peaks.plot);
	}

driver.n.peak.plots <- lapply(
	unique(driver.clusters$cluster),
	function(driver.cluster)
		create.n.peaks.plot(
			drivers = num.sig.peaks.df[num.sig.peaks.df$cluster == driver.cluster, 'driver'],
			num.sig.peaks.df = num.sig.peaks.df,
			xlabels.num.peaks,
			x.at.num.peaks,
			xlimits = c(0, max(num.sig.peaks.df$num.peaks.log) + .5),
			print.xlabels = driver.cluster == max(unique(driver.clusters$cluster))
			)
	);

plot.legend <- legend.grob(
	legends = c(mutation.legend, list(legend = m6A.enzyme.legend), list(legend = beta.legend)),
	layout = c(3, 1),
	label.cex = 0.55,
	between.col = 1,
	title.cex = 0.75,
	title.just = 'left',
	title.fontface = 'bold',
	size = 1.5
	);

#determine relative spacing of plots
driver.clust.sizes <- table(driver.clusters$cluster);
driver.clust.sizes[length(driver.clust.sizes)] <- driver.clust.sizes[length(driver.clust.sizes)] * 1.5;
y.heights <- determine.spacing(log2(driver.clust.sizes), available.space = 1, min.space = 0.1);
mutation.plot.width <- 0.1175;
m6A.enzyme.plot.width <- 0.044;
num.peaks.plot.width <- 0.1;
x.widths <- c(
	mutation.plot.width,
	m6A.enzyme.plot.width,
	determine.spacing(
			table(peak.clusters$cluster),
			available.space = 1 - sum(c(mutation.plot.width, m6A.enzyme.plot.width, num.peaks.plot.width)),
			min.space = 0.04
			),
		num.peaks.plot.width
		);

#combine plots
combined.plots <- unlist(lapply(
	unique(driver.clusters$cluster),
	function(driver.cluster)
		c(
			list(mutation.plots[[driver.cluster]]),
			list(m6A.enzyme.plots[[driver.cluster]]),
			unlist(lapply(cluster.heatmaps[[driver.cluster]], function(this.plot) this.plot), recursive = FALSE),
			list(driver.n.peak.plots[[driver.cluster]])
			)
	),
	recursive = FALSE
	);

#plotting parameters
plot.height <- 8;
y.spacing <- -0.51;
x.spacing <- -0.23;
plot.width <- 8;

#final composite plot
create.multipanelplot(
	combined.plots,
	filename = plot.path,
	layout.width = length(x.widths),
	layout.height = length(y.heights),
	size.units = 'in',
	resolution = 300,
	style = 'Nature',
	plot.objects.heights = y.heights,
	plot.objects.widths = x.widths,
	ylab.axis.padding = c(rep(0, length(x.widths) - 1), 1),
	# the last element of this vector creates some needed space between the plots and x-axis labels
	xlab.axis.padding = c(rep(0, length(y.heights) - 1), 2),
	y.spacing = y.spacing,
	x.spacing = x.spacing,
	layout.skip = rep(FALSE, length(y.heights) * length(x.widths)),
	width = plot.width,
	height = plot.height,
	bottom.padding = 0,
	top.padding = -1.5,
	top.legend.padding = 0,
	right.padding = 0,
	legend = list(
		bottom = list(
			fun = plot.legend
			)
		),
	);

### WRITE SESSION PROFILE TO FILE ######################################
if (!dir.exists(session.info.dir)) {
	dir.create(session.info.dir);
	}
save.session.profile(paste0(session.info.dir, timestamp, script.name, '.txt'));

#save config parameter information
write.table(
	t(data.frame(config)),
	paste0(session.info.dir, timestamp, script.name, '_config_params.txt'),
	sep = '\t',
	row.names = TRUE,
	col.names = FALSE,
	quote = FALSE
	);
