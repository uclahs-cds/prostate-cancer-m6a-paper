library(GenomicRanges);
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);

config <- 'D_E.config';

gene_id <- config$gene
gtf <- read_gtf(config$gtf, zero_based = FALSE);

# Gene Coordinates
genes.coords <- gtf[grep(gene_id, gtf$gene_id),
    c('chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'type')];
genes.coords <- data.frame(genes.coords);
rm(gtf);

m6a.seq.peaks <- read.delim('MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.peaks.tsv', row.names = 1)

m6a.seq <- read.delim('MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.Mask.tsv')

make.gene.peaks <- function(peaks) {
    if (nrow(peaks) > 0) {
        peaks <- data.frame(intersect(makeGRangesFromDataFrame(peaks), makeGRangesFromDataFrame(genes.coords[genes.coords$type == 'exon', ])))
        return(cbind(peaks, col = 'darkblue'))
        } else {return(NULL)}
    }
gene.peaks <- apply(m6a.seq[grep(gene_id, row.names(m6a.seq)), ], 2, 
    function(peak) make.gene.peaks(m6a.seq.peaks[grep(gene_id, row.names(m6a.seq.peaks))[as.logical(peak)], 1:6])
    )

sac.seq.gene <- read.delim(paste0(gene_id, '.sites.txt'))

overlaps.data <- list()
for (sample in grep('CPCG', names(m6a.seq), value = TRUE)){
    peaks <- grep(gene_id, row.names(m6a.seq))
    peaks <- peaks[as.logical(m6a.seq[peaks, sample])]
    overlaps.data[[sample]] <- data.frame(
        overlaps = length(which(overlapsAny(makeGRangesFromDataFrame(sac.seq.gene[sac.seq.gene$sample == sample, ]), makeGRangesFromDataFrame(m6a.seq.peaks[peaks, 1:6])))),
        sac.seq = length(makeGRangesFromDataFrame(sac.seq.gene[sac.seq.gene$sample == sample, ])),
        m6a.seq =  length(makeGRangesFromDataFrame(m6a.seq.peaks[peaks, 1:6]))
        )
    }
overlaps.data <- do.call(rbind, overlaps.data)

load(config$peak_coverage);

IP.coverage <- t(IP.coverage);
input.coverage <- t(input.coverage);

SAC.plot <- function(sample, gene.start = min(genes.coords$start), gene.end = max(genes.coords$end), xlimits = NULL) {
    sample.data <- sac.seq.gene[sac.seq.gene$sample == sample, ]
    rect.range = -10
    sample.data$bottom = 0.7*rect.range
    sample.data$top = 0.3*rect.range
    sample.data$start = sample.data$start - gene.start + 1
    sample.data$end = sample.data$end - gene.start + 1
    sample.data$col <- "darkblue"
    rect.data = data.frame(matrix(0, ncol = 4, nrow = (nrow(sample.data) * 2) - 1))
    colnames(rect.data) <- c('start', 'end', 'bottom', 'top')
    for(peak in 1:(nrow(sample.data))) {
        new.peak <- peak + (1 * (peak - 1))
        rect.data[new.peak, ] <- as.numeric(sample.data[peak, c('start', 'end', 'bottom', 'top')])
        if (new.peak < nrow(sample.data)) {
            rect.data[new.peak + 1, 'start'] <- as.numeric(sample.data[peak, 'end'] + 1)
        }
    }
    rect.data$start <- rect.data$start - 19
    rect.data$end <- rect.data$end + 19
    rect.data$col <- "darkblue"

    #only plot gene info
    gene.data = data.frame(x = 1:(gene.end - gene.start), min = 0, max = 100)
    if (is.null(xlimits)) {
        xlimits = c(0, max(gene.data$x)
        }
    sac.plt <- create.polygonplot(
        formula = NA ~ x,
        data = gene.data,
        max = gene.data$max,
        min = gene.data$min,
        # adjust axes limits
        xlimits = xlimits,
        ylimits = c(-8,-2),
        ylab.label = 'SAC-seq\nsites',
        ylab.cex = .5,
        xlab.cex = 0,
        yaxis.cex = 0,
        xaxis.cex = 0,
        yaxis.tck = 0,
        xaxis.tck = 0,
        main.cex = .7,
        main.y = .3,
        main.x = .6,
        # add fill colour
        col = 'seagreen3',
        border.col = 'white',
        xyline.lty = c(1, 0),
        add.border = TRUE,
        alpha = 0.2,
        # add middle line
        add.median = F,
        #median.lty = rep(c(3, 1), length(unique(pg.data$groups))/2),
        median.lwd = 1,
        #median.col = unlist(lapply(colours, rep, 2)),
        # Add rectangles
        add.rectangle = TRUE, #don't plot gene annotation for this peak region - too small
        xleft.rectangle = rect.data$start,
        xright.rectangle = rect.data$end,
        ybottom.rectangle = rect.data$bottom,
        ytop.rectangle = rect.data$top,
        col.rectangle = 'darkblue',
        alpha.rectangle = 1,
        top.padding = -0.5
        );
    return(sac.plt)
    }

samples <- grep('CPCG', names(m6a.seq), value = TRUE)
names(samples) = samples
coverage.list <- lapply(samples, function(sample)
    list(
        IP = IP.coverage[sample, , drop = FALSE],
        input = input.coverage[sample, , drop = FALSE]
        )
    )

process.coverage <- function(coverage.data, log.data = FALSE, input.IP.ratio = FALSE) {
    if (log.data) {
        coverage.data$IP <- log10(coverage.data$IP + 1);
        coverage.data$input <- log10(coverage.data$input + 1);
        }
    if (input.IP.ratio) {
        pg.data <- data.frame(
            'x' = 1:ncol(coverage.data$IP),
            'max' = c(apply(coverage.data$IP, 2, max) + 1) / c(apply(coverage.data$input, 2, max) + 1),
            'median' = c(apply(coverage.data$IP, 2, median) + 1) / c(apply(coverage.data$input, 2, median) + 1),
            'min' = rep(0, ncol(coverage.data$IP))
            );
        } else {
        pg.data <- data.frame(
            'x' = rep(1:ncol(coverage.data$IP), 2),
            'max' = c(apply(coverage.data$IP, 2, max), apply(coverage.data$input, 2, max)),
            'median' = c(apply(coverage.data$IP, 2, median), apply(coverage.data$input, 2, median)),
            'min' = rep(0, ncol(coverage.data$IP)*2),
            'groups' = c(
                rep('IP', ncol(coverage.data$IP)),
                rep('input', ncol(coverage.data$input))
                )
            );
        }
    return(pg.data)
    }

pg.data.conditions <- lapply(
    coverage.list,
    process.coverage,
    log.data = FALSE,
    input.IP.ratio = FALSE
    );

plot_gene <- function(pg.data, peak.region, gene.start = genes.coords$start[1],
    gene.end = max(genes.coords$end), main, max.data = max(pg.data$max) * 1.1,
    log.data = FALSE, yaxis = NULL, add.rect = T,
    xlimits = c(1, ifelse(is.null(pg.data$groups), nrow(pg.data), table(pg.data$groups)[[1]])),
    ylimits = NULL
    ) {

    # Gene Annotations
    if (add.rect) {
        rect.range <- -0.1*max.data
        peak.region$bottom = 0.7*rect.range
        peak.region$top = 0.3*rect.range
        peak.region$start = peak.region$start - gene.start + 1
        peak.region$end = peak.region$end - gene.start + 1
        rect.data = data.frame(matrix(0, ncol = 5, nrow = (nrow(peak.region) * 2) - 1))
        colnames(rect.data) <- c('start', 'end', 'bottom', 'top', 'col')
        for(peak in 1:(nrow(peak.region))) {
            new.peak <- peak + (1 * (peak - 1))
            rect.data[new.peak, c('start', 'end', 'bottom', 'top')] <- as.numeric(peak.region[peak, c('start', 'end', 'bottom', 'top')])
            rect.data[new.peak, 'col'] <- peak.region[peak, 'col']
            rect.data[new.peak + 1, 'start'] <- as.numeric(peak.region[peak, 'end'] + 1)
            rect.data[new.peak + 1, 'col'] <- 'white'
            if (new.peak < nrow(peak.region)) {
                rect.data[new.peak + 1, 'end'] <- as.numeric(peak.region[peak + 1, 'start'] - 1)
                } else {
                rect.data[new.peak + 1, 'end'] <- gene.end - gene.start
                }
        }
    } else {
        rect.data = data.frame(matrix(0, ncol = 5, nrow = 1))
        colnames(rect.data) <- c('start', 'end', 'bottom', 'top', 'col')
        rect.range <- 0
        rect.data$col <- "darkblue"
    }

    if (is.null(pg.data$groups)) {
        col <- 'seagreen3';
        median <- pg.data$max;
        } else {
        col <- c('darkorchid4', 'seagreen3');
        median <- pg.data$median;
        }
    if (is.null(ylimits)) {
        ylimits = c(rect.range, max(pg.data$max) * 1.1)
        }
    if (log.data) {
        yaxis <- auto.axis(10^max(ylimits), log.scaled=log.data, num.labels=3);
        } else {
        yaxis <- auto.axis(max(ylimits), log.scaled=log.data, num.labels=3);
        }
    # create the simple plot
    pg.plt <- create.polygonplot(
        formula = NA ~ x,
        #filename = paste0(path, gene_id, '_cell_line_peak_coverage.png'),
        data = pg.data,
        min = pg.data$min,
        max = pg.data$max,
        groups = pg.data$groups,
        # adjust axes limits
        xlimits = xlimits,
        ylimits = ylimits,
        ylab.cex = 1.5,
        xlab.cex = 0,
        yaxis.cex = .8,
        xaxis.cex = 0,
        yaxis.tck = 0,
        xaxis.tck = 0,
        yat = yaxis$at,
        yaxis.lab = yaxis$axis.lab,
        ylab.label = '',
        main = main,
        main.cex = .7,
        main.y = -.8,
        main.x = .6,
        # add fill colour
        col = col,
        border.col = 'white',
        xyline.lty = c(1, 0),
        add.border = TRUE,
        alpha = 0.2,
        # add middle line
        add.median = T,
        median = median,
        #median.lty = rep(c(3, 1), length(unique(pg.data$groups))/2),
        median.lwd = 1,
        #median.col = unlist(lapply(colours, rep, 2)),
        # Add rectangles
        add.rectangle = add.rect, #don't plot gene annotation for this peak region - too small
        xleft.rectangle = rect.data$start,
        xright.rectangle = rect.data$end,
        ybottom.rectangle = rect.data$bottom,
        ytop.rectangle = rect.data$top,
        col.rectangle = rect.data$col,
        alpha.rectangle = 1,
        top.padding = -0.5
        );
    return(pg.plt);
    }

plot.gene.schema <- function(genes.coords, xlimits = c(0, max(gene.data$x))) {
    gene <- genes.coords$gene_name[1]
    rect.range = -10
    # Gene Annotations
    rect.data = genes.coords
    gene.start = rect.data$start[rect.data$type == "gene"]
    rect.data$start = rect.data$start - gene.start + 1
    rect.data$end = rect.data$end - gene.start + 1
    rect.data = rect.data[rect.data$type %in% c("UTR", "exon"), c('start', 'end', 'type')]
    rect.data$bottom = ifelse(rect.data$type == "UTR", 0.55*rect.range, 0.9*rect.range)
    rect.data$top = ifelse(rect.data$type == "UTR", 0.45*rect.range, 0.1*rect.range)
    #fill in any gaps in annotation with thin/UTR bar
    rect.data = rbind(rect.data, data.frame('start' = 1, 'end' = max(rect.data$end), 'type' = 'gene', 'bottom' =  0.55*rect.range, 'top' = 0.45*rect.range))
    rect.data$col <- "darkblue"
    #only plot gene info
    gene.data = data.frame(x = 1:(max(rect.data$end)), min = 0, max = 100)
    gene.plot <- create.polygonplot(
        formula = NA ~ x,
        data = gene.data,
        min = min(gene.data$min),
        max = max(gene.data$max),
        # adjust axes limits
        xlimits = xlimits,
        ylimits = c(rect.range + 10,-10),
        ylab.label = gene,
        ylab.cex = .5,
        xlab.cex = 0,
        yaxis.cex = 0,
        xaxis.cex = 0,
        yaxis.tck = 0,
        xaxis.tck = 0,
        main.cex = .8,
        main.y = -.3,
        main.x = .6,
        # add fill colour
        border.col = 'white',
        xyline.lty = c(1, 0),
        add.border = TRUE,
        alpha = 0.2,
        # add middle line
        add.median = F,
        #median.lty = rep(c(3, 1), length(unique(pg.data$groups))/2),
        median.lwd = 1,
        #median.col = unlist(lapply(colours, rep, 2)),
        # Add rectangles
        add.rectangle = TRUE, #don't plot gene annotation for this peak region - too small
        xleft.rectangle = rect.data$start,
        xright.rectangle = rect.data$end,
        ybottom.rectangle = rect.data$bottom,
        ytop.rectangle = rect.data$top,
        col.rectangle = rect.data$col,
        alpha.rectangle = 1,
        top.padding = -0.5
        );
    return(gene.plot)
    }

if (genes.coords$gene_name[1] == 'AR') {
    xlimits = c(67711402, max(genes.coords$end)) - min(genes.coords$start)
    } else {
    xlimits = c(0, max(gene.data$x)) - min(genes.coords$start)
    }

gene.plot <- plot.gene.schema(genes.coords, xlimits = xlimits)
samples <- unique(sac.seq.gene$sample)[samples = unique(sac.seq.gene$sample) %in% names(gene.peaks)]
for (sample in samples) {

    SAC.plt <- SAC.plot(sample = sample,  xlimits = xlimits)
    m6A.plot <- plot_gene(
        main = sample,
        pg.data = pg.data.conditions[[sample]],
        peak.region = gene.peaks[[sample]],
        gene.start = genes.coords$start[1],
        add.rect = !is.null(gene.peaks[[sample]]),
        xlimits = xlimits
        )

    peak.plots <- list(m6A.plot, SAC.plt, gene.plot)
    create.multipanelplot(
        plot.objects = peak.plots,
        filename = paste0(genes.coords$gene_name[1], '_', sample, '_gene_level.pdf'),
        layout.height = length(peak.plots),
        plot.objects.heights = c(.7, .15, .15),
        layout.width = 1,
        y.spacing = -.8,
        ylab.axis.padding = -0.25,
        height = 3,
        width = 4,
        top.padding = -2.5,
        bottom.padding = -2.5,
        );
    }
