library(GenomicRanges)
library(BoutrosLab.plotting.general)
library(VennDiagram)

#load tumor and normal BED m6A peak files
normal.bed <- Sys.glob('input_data/M6A_NormalTumour/*/peak.bed')
tumor.bed <- Sys.glob('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/peak_fitting/F_MeTPeak_Peaks_SS/*/peak.bed')
names(tumor.bed) <- basename(dirname(tumor.bed))
tumor.joint.peaks <- read.delim('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/processed_data/GRCh38/Z_Matrices/PeakCounts_SoftMask/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakCounts_SoftMask.Adjusted.tsv')
#keep only these samples
tumor.samples <- names(tumor.joint.peaks)
tumor.bed <- tumor.bed[tumor.samples]

load.bed.files <- function(bed.files) {
    names(bed.files) <- basename(dirname(bed.files))
    beds <- lapply(bed.files, read.delim, check.names = FALSE)
    #filter by peak number
    beds <- beds[which(lapply(beds, nrow) > 2000)]
    GRs <- GRangesList(lapply(beds, makeGRangesFromDataFrame, seqnames.field = '# chr', start.field = 'chromStart', end.field = 'chromEnd'))
    return(GRs)
    }

compute.intersections <- function(GR.x, GR.y, return.diag) {
    intersections <- sapply(GR.x, function(y) sapply(GR.y, function(x) length(which(overlapsAny(x, y))) / length(x)))
    intersections <- round(intersections, 2);
    intersections[lower.tri(intersections, diag = return.diag)]
    }

normal.GRs <- load.bed.files(normal.bed)
tumor.GRs <- load.bed.files(tumor.bed)
norm.norm <- compute.intersections(normal.GRs, normal.GRs, return.diag = FALSE)
tumor.tumor <- compute.intersections(tumor.GRs, tumor.GRs, return.diag = FALSE)
norm.tumor <- compute.intersections(normal.GRs, tumor.GRs, return.diag = TRUE)

to.plot <- data.frame(
    proportion = c(norm.norm, norm.tumor, tumor.tumor),
    comparison = c(rep('Normal-Normal', length(norm.norm)), rep('Normal-Tumor', length(norm.tumor)), rep('Tumor-Tumor', length(tumor.tumor)))
    )

create.boxplot(
    comparison ~ proportion,
    data = to.plot,
    filename = 'results/normal_samples/peak_overlap_boxplot.pdf',
    height = 2.5,
    width = 3,
    add.stripplot = TRUE,
    xaxis.fontface = 'bold',
    yaxis.fontface = 'bold',
    ylab.label = '',
    xlab.label = 'Proportion Shared Peaks',
    points.pch = 21,
    main.cex = .75,
    xlab.cex = .6,
    ylab.cex = .6,
    xaxis.cex = .5,
    yaxis.cex = .5,
    left.padding = 2,
    xaxis.tck = 0.2,
    yaxis.tck = 0.3
    );

venn.diagram(
    x = list(tumor.peaks, normal.peaks),
    category.names = c('Tumor', 'Normal'),
    filename = 'results/normal_samples/normal_tumor_peak_venn.png',
    output = TRUE,
    hyper.test = FALSE,
    lower.tail = FALSE, #test the probability of obtaining at least the observed intersection,
    total.population = length(c(normal.peaks, tumor.peaks)),
    # Output features
    imagetype = 'png',
    resolution = 300,
    # Circles
    lwd = .5,
    #lty = 'blank',
    fill = c('darkorchid4', 'darkgreen'),
    # Numbers
    cex = .9,
    fontface = 'bold',
    fontfamily = 'ariel',
    # Set names
    cat.cex = 0.9,
    cat.fontface = 'bold',
    cat.default.pos = 'outer',
    cat.fontfamily = 'ariel',
    cat.pos = c(200, 160),
    #margin = 3,
    ext.dist = c(-0.1, 0.02),
    #ext.length = rep(.9, 2),
    #sub.pos = c(0.5, 0.55),
    main.fontfamily = 'ariel',
    units = 'in',
    height = 4.5,
    width = 4.5
    );
