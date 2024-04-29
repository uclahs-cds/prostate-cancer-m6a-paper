library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
source('functions/diverging_heatmap.R');
source('functions/FDR_barplots.R');

script.name <- 'survival_peak_forestplot.R';
timestamp.now <- format(Sys.time(), '%Y%m%d-%H%M%S');
session.info.dir <- './session_info/';
config.path <- '/hot/users/rhughwhite/projects/PRAD-000080-N6Methyl/code/project-ProstateCancer-m6A/survival/config/survival_peak_forestplot.config';
config <- read.config.file(config.path);

setwd(config$working.dir);
output.dir <- config$output.dir;

survival.results <- read.delim(config$survival.results);
survival.results <- survival.results[which(survival.results$qvalue < .1), ];
survival.results <- survival.results[order(survival.results$hr, decreasing = FALSE), ];


survival.results$log.HR <- log2(survival.results$hr);
survival.results$log.HR.high <- log2(survival.results$upper.95);
survival.results$log.HR.low <- log2(survival.results$lower.95);

survival.results$label <- survival.results$gene.name;
x.at.values <- seq(0, 4, 2);
xlimits <- log2(c(0.95, 2^5));
xaxis.cex = .45;
yaxis.cex = .65;

p.value.barplot <- FDR.barplot(
    FDR = survival.results$pvalue,
    log10 = TRUE,
    FDR.threshold = 0.05,
    FDR.intervals = -log10(c(1, .1, .001, .00001, .000001)),
    first.xaxis.lab = '0',
    plot.xaxis.labels = TRUE,
    xlab.label = 'P',
    xlab.cex = .65,
    xaxis.cex = xaxis.cex
    );

x.at.labels <- sapply(x.at.values, function(power) as.expression(bquote(bold(.('2') ^.(power)))));
forest <- create.scatterplot(
    formula = 1:nrow(survival.results) ~ log.HR,
    data = survival.results,
    ylab.label = 'Peak',
    xlab.label = 'Hazard Ratio',
    yaxis.lab = survival.results$label,
    yat = 1:nrow(survival.results),
    ylimits = c(0.5, nrow(survival.results) + .5), #need to add .5 spacing to line up with barplot
    xlimits = xlimits,
    xaxis.lab = x.at.labels,
    xat = log2(sapply(x.at.values, function(x) 2^x)),
    xlab.cex = .65,
    xaxis.cex = .6,
    yaxis.cex = yaxis.cex,
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

plots <- list(gene.heatmap, forest, p.value.barplot);
plot.objects.widths = c(0.25, 0.45, 0.3);
width <- 3;
height <- 5;
x.spacing <- -.6;


filename <- paste0(output.dir, 'survival_peak_forestplot_all_peaks.pdf');
final.plot <- create.multipanelplot(
    plots,
    layout.width = length(plots),
    layout.height = 1,
    size.units = 'in',
    resolution = 300,
    style = 'Nature',
    plot.objects.widths = plot.objects.widths,
    x.spacing = x.spacing,
    bottom.padding = -2,
    top.padding = 0,
    left.padding = -2,
    top.legend.padding = -1,
    right.legend.padding = 0.5,
    right.padding = -.5
    );

pdf(filename, width = width, height = height);
print(final.plot);
dev.off();

### WRITE SESSION PROFILE TO FILE ######################################
if (!dir.exists(session.info.dir)) {
    dir.create(session.info.dir);
    }
save.session.profile(paste0(session.info.dir, timestamp.now, script.name, '.txt'));

#save config parameter information
write.table(
    t(data.frame(config)),
    paste0(session.info.dir, timestamp.now, script.name, '_config_params.txt'),
    sep = '\t',
    row.names = TRUE,
    col.names = FALSE,
    quote = FALSE
    );
