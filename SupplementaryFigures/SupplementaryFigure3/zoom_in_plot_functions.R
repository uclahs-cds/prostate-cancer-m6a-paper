library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

plot.my.scatterplot <- function(
	var1,
	var2,
	labels,
	file,
	main.cex = 0.7,
	alpha = 0.75,
	xlab.label = expression(bold('m')^bold('6') * bold('A abundance')),
	col = 'black',
	key.y = 1,
	height = 2.5,
	width = 2.5,
	...
	) {
	toplot <- data.frame(
		var1 = as.numeric(var1),
		var2 = as.numeric(var2),
		stringsAsFactors = FALSE
		);

	plot <- create.scatterplot(
		formula =  var1 ~ var2,
		data = toplot,
		height = 2.5,
		width = 2.5,
		ylab.label = labels[1],
		xlab.label = xlab.label,
		main = labels[2],
		main.cex = main.cex,
		col = col,
		ylab.cex = 0.65,
		xlab.cex = 0.65,
		xaxis.cex = 0.65,
		yaxis.cex = 0.65,
		xaxis.tck = c(0.5, 0),
		yaxis.tck = c(0.5, 0),
		cex = 0.3,
		alpha = alpha,
		type = c('p','r'),
		lwd = 1.5,
		lty = 2,
		resolution = 300,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = get.corr.key(
						x = toplot$var1,
						y = toplot$var2,
						label.items = c('spearman', 'spearman.p'),
						alpha.background = 0,
						key.cex = 0.5,
						)
					),
					x = 0.45,
					y = key.y
				)
			),
		...
		);
	plot$par.settings$axis.components <- list(bottom = list(pad1 = 1, pad2 = -0.5), left = list(pad1 = 1, pad2 = -0.5));
	png(file, height = height, width = width, units = 'in', res = 300);
	print(plot);
	dev.off();
	}

plot.my.boxplot <- function(var1, var2, labels, test, file, colour = 'black', status.one, status.two, main.cex = .7, ylab.label = expression(bold('m')^bold('6') * bold('A abundance')), max.value = NA, pch = 21) {

	toplot <- data.frame(
		var1 = as.factor(var1),
		var2 = as.numeric(var2),
		stringsAsFactors = FALSE
		);

	if (test == 'utest') {
		p.value <- wilcox.test(
			x = as.numeric(var2[which(var1 == status.one)]),
			y = as.numeric(var2[which(var1 == status.two)])
			)$p.value;
		x.pos <- 0.45;
		}

	if (test == 'aov') {
		aov.result <- summary(
			aov(
				formula = var2 ~ var1,
				data = toplot
				)
			);
		p.value <- aov.result[[1]][['Pr(>F)']][1];
		x.pos <- 0.45;
		}

	if (round(p.value, digits = 3) == 0) {
		your.number <- list(
			base = scientific.notation(p.value, digits = 2, type = 'list')$base,
			exponent = scientific.notation(p.value, type = 'list')$exponent
			);
		my.legend <- as.expression(substitute(P == paste(base %*% 10^exponent), your.number));
	} else {
		my.legend <- paste0('P = ', round(p.value, digits = 3));
		}

	if (is.na('max.value')) {
	} else {
		max.value <- round(max(var2), -2);
		}
	if (max.value == 0) {
		max.value <- max(var2);
		yat <- seq(0, max.value, 10);
	} else if (max.value > 0 & max.value < 100) {
		yat <- seq(0, max.value, 10);
	} else if (max.value > 100 & max.value <= 500) {
		yat <- seq(0, max.value, 100);
	} else if (max.value > 500 & max.value < 2000) {
		yat <- seq(0, max.value, 200);
	} else if (max.value == 15100) {
		yat <- seq(0, max.value, 5000);
	} else {
		yat <- seq(0, max.value, 1000);
		}

	create.boxplot(
		formula = var1 ~ var2,
		data = toplot,
		height = 2.64,
		width = 2.76,
		filename = file,
		add.stripplot = TRUE,
		ylab.label = labels[1],
		xlab.label = ylab.label,
		main = labels[2],
		xlimits = c(0, max.value),
		xat  = yat,
		points.col = colour,
		points.pch = pch,
		points.cex = 0.4,
		main.cex = main.cex,
		xaxis.tck = c(0.5, 0),
		yaxis.tck = c(0.5, 0),
		xaxis.cex = 0.65,
		yaxis.cex = 0.65,
		ylab.cex = 0.65,
		xlab.cex = 0.65,
		resolution = 300,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = list(
					text = list(
						lab = my.legend,
						cex = 0.5,
						fontface = 'plain'
						)
					)
				),
				x = x.pos,
				y = 0.1
				)
			)
		);
	}
