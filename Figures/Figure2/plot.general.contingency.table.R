library(BoutrosLab.plotting.general);
library(chisq.posthoc.test);

plot.general.contingency.table <- function(
	plot.data = plot.data,
	filename = filename,
	plot.title = '',
	xlabel = '',
	ylabel = '',
	left.padding = 1,
	width = 5,
	height = 5,
	panel.height = 0.9,
	panel.width = 0.9,
	xaxis.rot = 0,
	scheme = 'input',
	text.cex = 0.95,
	main.cex = 1,
	ylab.cex = 1,
	xlab.cex = 1,
	xaxis.cex = 0.95,
	yaxis.cex = 0.95,
	axes.lwd = 0.5,
	row.lwd = 0.5,
	col.lwd = 0.5,
	...
	) {

	# check if plot.data is a table
	if (class(plot.data) != 'table') {
		return('Check your input contingency table\n');
		}

	results <- chisq.test(plot.data);
	posthoc.results <- chisq.posthoc.test(
		x = plot.data,
		method = 'none'
		);
	pvalue.table <- posthoc.results[which(posthoc.results$Value == 'p values'), -match(c('Dimension', 'Value'), colnames(posthoc.results))];
	residual.table <- posthoc.results[which(posthoc.results$Value == 'Residuals'), -match(c('Dimension', 'Value'), colnames(posthoc.results))];

	colour.matrix <- matrix(
		data = 'white',
		nrow = nrow(pvalue.table),
		ncol = ncol(pvalue.table)
		);
	colour.matrix[which(pvalue.table < 0.05 & residual.table < 0)] <- 'forestgreen';
	colour.matrix[which(pvalue.table < 0.05 & residual.table > 0)] <- 'darkorchid';

	if (round(results$p.value, digits = 2) == 0) {
		your.number <- list(
			base = scientific.notation(results$p.value, digits = 2, type = 'list')$base,
			exponent = scientific.notation(results$p.value, type = 'list')$exponent
			);
		main.text <- substitute(bold(P == paste(base %*% 10^exponent)), your.number);
		main.text <- bquote(bold("P"==.(your.number$base)~"x"~.(as.character(10))^.(as.character((your.number$exponent)))));
	} else {
		main.text <- paste0('P = ', round(results$p.value, digits = 3));
		}

	if (scheme == 'count') {
		x <- plot.data / sum(plot.data);
		colour.scheme <- c('white', 'dodgerblue');
		input.colours <- FALSE;
		at <- seq(0, 1, 0.1);
		final.colour <- 'dodgerblue';
	} else {
		x <- colour.matrix;
		colour.scheme <- NULL;
		input.colours <- TRUE;
		at <- NULL;
		final.colour <- 'snow';
		}

	# create heatmap
	main.heatmap <- create.heatmap(
		x = x,
		clustering = 'none',
		input.colours = input.colours,
		at = at,
		colour.scheme = colour.scheme,
		cell.text = data.frame(plot.data)$Freq,
		text.cex = text.cex,
		text.fontface = 2,
		grid.row = TRUE,
		grid.col = TRUE,
		xlab.label = xlabel,
		ylab.label = ylabel,
		ylab.cex = ylab.cex,
		xlab.cex = xlab.cex,
		xat = 1:ncol(plot.data),
		yat = 1:nrow(plot.data),
		xaxis.lab = colnames(plot.data),
		yaxis.lab = rownames(plot.data),
		xaxis.cex = xaxis.cex,
		yaxis.cex = yaxis.cex,
		axes.lwd = axes.lwd,
		row.lwd = row.lwd,
		col.lwd = col.lwd,
		col.pos = rep(1:ncol(plot.data), each = nrow(plot.data)),
		row.pos = rep(nrow(plot.data):1, times = ncol(plot.data)),
		print.colour.key = FALSE,
		same.as.matrix = TRUE,
		xaxis.rot = xaxis.rot,
		use.legacy.settings = FALSE,
		...
		);

	# row total heatmap
	row.total <- create.heatmap(
		x = data.frame(rowSums(plot.data), rowSums(plot.data), rowSums(plot.data)) / sum(plot.data),
		clustering = 'none',
		colour.scheme = c('white', final.colour),
		colour.alpha = 1,
		at = at,
		cell.text = rev(rowSums(plot.data)),
		text.cex = text.cex,
		text.fontface = 2,
		xaxis.cex = 0,
		axes.lwd = axes.lwd,
		row.lwd = row.lwd,
		col.lwd = col.lwd,
		grid.row = TRUE,
		grid.col = FALSE,
		col.pos = rep(2, nrow(plot.data)),
		row.pos = 1:nrow(plot.data),
		print.colour.key = FALSE,
		same.as.matrix = TRUE,
		xaxis.tck = c(0, 0),
		use.legacy.settings = FALSE,
		...
		);

	# column total heatmap
	col.total <- create.heatmap(
		x = data.frame(colSums(plot.data)) / sum(plot.data),
		clustering = 'none',
		colour.scheme = c('white', final.colour),
		colour.alpha = 1,
		at = at,
		cell.text = colSums(plot.data),
		text.cex = text.cex,
		text.fontface = 2,
		main.cex = main.cex,
		yaxis.cex = 0,
		axes.lwd = axes.lwd,
		row.lwd = row.lwd,
		col.lwd = col.lwd,
		grid.row = FALSE,
		grid.col = TRUE,
		col.pos = 1:ncol(plot.data),
		row.pos = rep(1.5, ncol(plot.data)),
		print.colour.key = FALSE,
		yaxis.tck = c(0, 0),
		xaxis.top.tck = c(0, 0),
		same.as.matrix = TRUE,
		use.legacy.settings = FALSE,
		...
		);

	# create multiplot
	create.multipanelplot(
		plot.objects = list(main.heatmap, row.total, col.total),
		layout.width = 2,
		layout.height = 2,
		main = as.expression(main.text),
		main.cex = main.cex,
		filename = filename,
		resolution = 300,
		plot.objects.heights = c(panel.height, 0.2),
		plot.objects.widths = c(panel.width, 0.18),
		y.spacing = -0.5, 
		x.spacing = -0.5,
		layout.skip = c(FALSE, FALSE, FALSE, TRUE),
		xlab.axis.padding = 1,
		width = width,
		height = height
		);
	}
