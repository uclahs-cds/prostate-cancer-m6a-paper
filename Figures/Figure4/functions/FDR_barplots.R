### volcano_plots.R ############################################
# function to create volcano plots

library(BoutrosLab.plotting.general);

# FDR = vector of FDR values
# FDR.threshold = where to plot vertical line indicating FDR threshold of interest
# FDR.intervals = x axis values to label
# log10 = whether to transform FDR values
# plot.xaxis.labels = whether to plot x axis label and tick labels
# xlab.label = xlab.label
FDR.barplot <- function(
    FDR,
    FDR.threshold = 0.05,
    FDR.min = NULL,
    log10 = FALSE,
    plot.xaxis.labels = TRUE,
    xlab.label = 'FDR',
    xaxis.tck = c(.3, 0),
    yaxis.tck = 0,
    xaxis.cex = .4,
    xlab.cex = .7,
    FDR.intervals = NULL,
    first.xaxis.lab = NULL,
    FDR.seq.by = 2,
    ...
    ) {

    if (log10) {
        FDR <- data.frame(
            x = 1:length(FDR),
            FDR = -log10(FDR)
            );
        FDR.threshold <- -log10(FDR.threshold);
        if (is.null(FDR.min)) {
            FDR.min <- max(FDR$FDR + .5, na.rm = TRUE);
            }
        if (is.null(FDR.intervals)) {
            FDR.intervals <- seq(0, FDR.min, by = FDR.seq.by);
            }
        FDR.notation <- sapply(FDR.intervals, function(power) as.expression(bquote(bold(.('10') ^.(paste0('-', power))))));
        xlimits <- c(-0.05, FDR.min);
        } else {
        FDR <- data.frame(
            x = 1:length(FDR),
            FDR = 1 - FDR
            );
        FDR.threshold <- 1 - FDR.threshold;
        if (is.null(FDR.intervals)) {
            FDR.intervals = c(0, .5, .95);
            }
        FDR.intervals <- 1 - FDR.intervals;
        FDR.notation <- as.character(1 - FDR.intervals);
        xlimits <- c(-0.01, 1)
        }
    if (!is.null(first.xaxis.lab)) {
        FDR.notation[1] <- as.expression(bquote(bold(.(first.xaxis.lab))));
        }
    if (!plot.xaxis.labels) {
        FDR.intervals <- '';
        xlab.label <- '';
        FDR.notation <- '';
        }

    FDR.barplot <- create.barplot(
        as.formula('x ~ FDR'),
        data = FDR,
        plot.horizontal = TRUE,
        abline.v = FDR.threshold,
        abline.col = 'grey',
        border.lwd = 1,
        ylab.axis.padding = 0,
        left.padding = 0,
        xaxis.cex = xaxis.cex,
        ylab.label = '',
        yaxis.cex = 0,
        yaxis.lab = '',
        xat = FDR.intervals,
        xlimits = xlimits,
        xlab.label = bquote(bold(.(xlab.label))),
        xlab.cex = xlab.cex,
        xaxis.tck = xaxis.tck,
        yaxis.tck = yaxis.tck,
        xaxis.lab = FDR.notation,
        ...
        );
    return(FDR.barplot);
    }