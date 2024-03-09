library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

source('functions/diverging_heatmap.R');
source('functions/standard_barplot.R');
source('functions/discrete_colours_heatmap.R');
source('functions/CNA_heatmap.R');
source('functions/NA_legend.R');
source('functions/FDR_barplots.R');

script.name <- 'hypoxia_peak_heatmap.R';
timestamp.now <- format(Sys.time(), '%Y%m%d-%H%M%S');
session.info.dir <- './session_info/';
config.path <- '/hot/user/rhughwhite/projects/PRAD-000080-N6Methyl/code/project-ProstateCancer-m6A/hypoxia/configs/hypoxia_peak_heatmap.config';
config <- read.config.file(config.path);

setwd(config$working.dir);
output.dir <- config$output.path;

#generate heatmap of m6A peaks correlated with hypoxia
hyp.m6A <- read.delim(config$hypoxia.m6A);
hyp.m6A <- hyp.m6A[hyp.m6A$qvalue < config$FDR.threshold, ];
#load peak file
m6A.peaks <- read.delim(config$m6A.peaks);

if (config$peak.order.method == 'rho') {
    hyp.m6A <- hyp.m6A[order(hyp.m6A[config$peak.order.method], decreasing = TRUE), ];
    m6A.peaks <- m6A.peaks[hyp.m6A$peak.id, ];
    } else if (config$peak.order.method == 'clustering') {
    m6A.peaks <- m6A.peaks[hyp.m6A$peak.id, ];
    m6A.peaks.dend <- create.dendrogram(t(m6A.peaks)); #cluster with correlation distance and diana method
    peak.clusters <- cutree(as.hclust(m6A.peaks.dend), k = 2)[order.dendrogram(m6A.peaks.dend)];
    peak.clusters <- lapply(unique(peak.clusters), function(cluster) which(peak.clusters == cluster));
    m6A.peaks <- m6A.peaks[order.dendrogram(m6A.peaks.dend), ];
    hyp.m6A <- hyp.m6A[match(row.names(m6A.peaks), hyp.m6A$peak.id), ];
    }

buffa <- read.delim(
    file = config$hypoxia,
    as.is = TRUE
    );
#load omics data list - used in mutual information analysis
load(config$omics.data);

RNA <- read.delim(config$RNA); #using array data for RNA abundance to increase sample overlap with m6A

#only use samples with data in at least one other data type
buffa <- buffa[buffa$patient.IDs %in% names(m6A.peaks) & buffa$patient.IDs %in% names(RNA), ];
buffa <- buffa[order(buffa$scores, decreasing = FALSE), ];
#match columns with hypoxia, introduce NA if no matching sample
m6A.peaks <- t(t(m6A.peaks)[
    match(buffa$patient.IDs, colnames(m6A.peaks)),
    ]);
#list of m6A enzymes/RBPs
m6A.regulators <- read.delim('/hot/user/rhughwhite/projects/PRAD-000080-N6Methyl/input_data/2022-09-22_m6A.Enzymes_correct_IGF2BP.tsv');
m6A.regulators <- m6A.regulators$GeneName[m6A.regulators$ShortenedList == 1];

#m6A enzymes vs hypoxia
m6A.enzyme.CNA <- read.delim('/hot/users/rhughwhite/projects/PRAD-000080-N6Methyl/results/hypoxia/hypoxia_CNA_anova.txt');
m6A.enzyme.CNA <- m6A.enzyme.CNA[m6A.enzyme.CNA$gene %in% m6A.regulators, ];
if (as.logical(config$use.qvalue)) {
    m6A.enzyme.CNA$qvalue <- p.adjust(m6A.enzyme.CNA$pvalue, method = 'fdr');
    } else {
    m6A.enzyme.CNA$qvalue <- m6A.enzyme.CNA$pvalue;
    }
m6A.enzyme.CNA <- m6A.enzyme.CNA[order(m6A.enzyme.CNA$qvalue), ];

#RNA-hypoxia correlation
RNA.hypoxia <- read.delim('/hot/user/rhughwhite/projects/PRAD-000080-N6Methyl/results/hypoxia/hypoxia_RNA_correlations.txt');
m6A.enzyme.RNA <- RNA.hypoxia[RNA.hypoxia$gene %in% m6A.regulators, ];
if (!as.logical(config$use.qvalue)) {
    m6A.enzyme.RNA$qvalue <- m6A.enzyme.RNA$pvalue;
    }
m6A.enzyme.RNA <- m6A.enzyme.RNA[order(m6A.enzyme.RNA$qvalue), ];

#Format RNA data
RNA <- RNA[RNA$Symbol_UCSC %in% m6A.enzyme.RNA$gene, ];
#match columns with hypoxia, introduce NA if no matching sample
RNA <- t(t(RNA)[
    match(buffa$patient.IDs, colnames(RNA)),
    match(m6A.enzyme.RNA$gene, RNA$Symbol_UCSC)
    ]);
RNA <- apply(RNA, 2, as.numeric);
row.names(RNA) <- m6A.enzyme.RNA$gene;
#match columns with hypoxia, introduce NA if no matching sample
omics.data$CNA <- t(t(omics.data$CNA)[
    match(buffa$patient.IDs, colnames(omics.data$CNA)),
    m6A.enzyme.CNA$gene
    ]);

#hypoxia peaks pathway enrichment
gene.cluster.pathways <- read.delim('/hot/users/rhughwhite/projects/PRAD-000080-N6Methyl/results/hypoxia/hypoxia_m6A_pathway_genes.txt', check.names = FALSE);
pathway.annotation <- read.delim('/hot/users/rhughwhite/projects/PRAD-000080-N6Methyl/results/hypoxia/hypoxia_m6A_pathway_enrichment.txt', check.names = FALSE);

names(gene.cluster.pathways) <- pathway.annotation$term_name[match(names(gene.cluster.pathways), pathway.annotation$term_id)];
#shorten term for plot
names(gene.cluster.pathways) <- sub(
    'regulation of transcription by RNA polymerase II',
    'regulation of transcription by RNA pol. II',
    names(gene.cluster.pathways)
    );
#cluster peak data, do outside of heatmap function so that ordering can be used for other plots

RNA.hypoxia <- RNA.hypoxia[match(hyp.m6A$gene.name, RNA.hypoxia$gene), ];
gene.cluster.pathways <- gene.cluster.pathways[match(hyp.m6A$gene.name, row.names(gene.cluster.pathways)), ];

peak.gene.annotation <- data.frame(
    'RNA-hypoxia correlation' = RNA.hypoxia$qvalue < config$FDR.threshold,
    gene.cluster.pathways,
    check.names = FALSE,
    row.names = hyp.m6A$peak.id
    );
#barplot of hypoxia score
hypoxia.barplot <- create.standard.barplot(
    data = buffa,
    x = 'scores',
    ylab.label = 'Hypoxia Score',
    ylab.cex = .6,
    yaxis.cex = .4,
    yaxis.tck = c(.3, 0)
    );
hypoxia.barplot$par.settings$axis.components <- list(bottom = list(pad1 = .5), left = list(pad1 = .5));

m6A.heatmap.legends <- lapply(peak.clusters, function(cluster) 
    create.diverging.heatmap(
        data = m6A.peaks[names(cluster), ],
        same.as.matrix = TRUE,
        scale.data = 'row',
        resolution = 300,
        cluster.dimensions = 'none',
        plot.dendrograms = 'none',
        key.title = expression(bold(underline('Z'[m^6*A~abundance]))),
        yaxis.tck = c(0, 0),
        key.min = -2,
        key.max = 2
        ));

peak.gene.annotation.heatmap.legends <- lapply(peak.clusters, function(cluster)
    create.discrete.colour.heatmap(
        x = peak.gene.annotation[names(cluster), ],
        colour.scheme = as.list(setNames(c('cadetblue2', 'darkorange1', 'olivedrab2', 'black', 'tomato2'), names(peak.gene.annotation))),
        cluster.rows = FALSE,
        key.title = 'Hypoxia correlated peaks',
        same.as.matrix = TRUE
        ));

CNA.heatmap.legend <- create.CNA.heatmap(
    CNA.data = omics.data$CNA,
    cluster.dimensions = 'none',
    plot.dendrograms = 'none',
    yaxis.lab = row.names(omics.data$CNA),
    yaxis.cex = .35,
    yaxis.tck = c(.3, 0)
    );

CNA.FDR.barplot <- FDR.barplot(
    FDR = rev(m6A.enzyme.CNA$qvalue),
    FDR.threshold = 0.1,
    log10 = TRUE,
    #FDR.intevals = c(1, .5, .1),
    plot.xaxis.labels = TRUE,
    xlab.label = ifelse(as.logical(config$use.qvalue), 'Q', 'P'),
    xlab.cex = .5,
    xaxis.cex = .55
    );

RNA.heatmap.legend <- create.diverging.heatmap(
    data = RNA,
    same.as.matrix = TRUE,
    scale.data = 'row',
    cluster.dimensions = 'none',
    plot.dendrograms = 'none',
    key.title = expression(bold(underline('Z'[RNA~abundance]))),
    colours = c('blue', 'white', 'red'),
    yaxis.lab = row.names(RNA),
    yaxis.cex = .35,
    yaxis.tck = c(.3, 0)
    );
RNA.heatmap.legend$heatmap$par.settings$axis.components <- list(bottom = list(pad1 = .5), left = list(pad1 = .5));

RNA.FDR.barplot <- FDR.barplot(
    FDR = rev(m6A.enzyme.RNA$qvalue),
    FDR.threshold = 0.1,
    log10 = TRUE,
    plot.xaxis.labels = FALSE
    );

plot.legend <- legend.grob(
    legends = list(
        legend = m6A.heatmap.legends[[1]]$legend,
        legend = peak.gene.annotation.heatmap.legends[[1]]$legend,
        legend = RNA.heatmap.legend$legend,
        legend = CNA.heatmap.legend$legend
        ),
    label.cex = 0.5,
    #title.cex = 0.6,
    title.just = 'left',
    title.fontface = 'bold',
    size = 1.6
    );

plots <- list(
    hypoxia.barplot,
    m6A.heatmap.legends[[1]]$heatmap,
    peak.gene.annotation.heatmap.legends[[1]]$heatmap,
    m6A.heatmap.legends[[2]]$heatmap,
    peak.gene.annotation.heatmap.legends[[2]]$heatmap,
    RNA.heatmap.legend$heatmap,
    RNA.FDR.barplot,
    CNA.heatmap.legend$heatmap,
    CNA.FDR.barplot
    );

if(config$FDR.threshold == 0.1) {
    height <- 8;
    plot.objects.heights <- c(1.2, 6.25,
        sapply(peak.clusters, function(cluster) (length(cluster) / length(unlist(peak.clusters)) * 2)),
        2.25);
    } else {
    height <- 6;
    plot.objects.heights <- c(1.2,
        sapply(peak.clusters, function(cluster) (length(cluster) / length(unlist(peak.clusters)) * 2)),
        2.25, 2.25);
    }

create.multipanelplot(
    plots,
    layout.skip = c(FALSE, TRUE, rep(FALSE, length(plots) -1)),
    filename = paste0('results/hypoxia/hypoxia_m6A_peak_summary_FDR_', config$FDR.threshold, 'peak_order_', config$peak.order.method, '.pdf'),
    layout.width = 2,
    layout.height = 5,
    size.units = 'in',
    resolution = 300,
    style = 'Nature',
    plot.objects.heights = plot.objects.heights,
    plot.objects.widths = c(5, 1),
    #ylab.axis.padding = c(rep(0, length(x.widths) - 1), 1),
    # the last element of this vector creates some needed space between the plots and x-axis labels
    #xlab.axis.padding = c(rep(0, length(y.heights) - 1), 2),
    y.spacing = c(-3, -5, -5),
    x.spacing = -2.5,
    ylab.axis.padding = c(-3, 0),
    width = 6.75,
    height = height,
    bottom.padding = -2,
    top.padding = 0,
    left.padding = -1.75,
    top.legend.padding = -1,
    right.legend.padding = 0.5,
    right.padding = 1,
    legend = list(
        right = list(
            fun = plot.legend
            )
        ),
    );

RNA.FDR.barplot <- FDR.barplot(
    FDR = rev(m6A.enzyme.RNA$qvalue),
    log10 = TRUE,
    FDR.threshold = 0.1,
    FDR.intervals = -log10(c(.1, 10^-5, 10^-10, 10^-15)),
    FDR.seq.by = 5,
    plot.xaxis.labels = TRUE,
    xlab.label = ifelse(as.logical(config$use.qvalue), 'Q', 'P'),
    xlab.cex = .5,
    xaxis.cex = .35
    );

RNA.FDR.barplot$par.settings$axis.components <- list(bottom = list(pad1 = .5), left = list(pad1 = .5));

plots <- list(
    hypoxia.barplot,
    m6A.heatmap.legends[[1]]$heatmap,
    peak.gene.annotation.heatmap.legends[[1]]$heatmap,
    m6A.heatmap.legends[[2]]$heatmap,
    peak.gene.annotation.heatmap.legends[[2]]$heatmap,
    RNA.heatmap.legend$heatmap,
    RNA.FDR.barplot
    );

plot.legend <- legend.grob(
    legends = list(
        legend = m6A.heatmap.legends[[1]]$legend,
        legend = peak.gene.annotation.heatmap.legends[[1]]$legend,
        legend = RNA.heatmap.legend$legend
        ),
    label.cex = 0.5,
    title.cex = .75,
    title.just = 'left',
    title.fontface = 'bold',
    size = 1.6
    );

if(config$FDR.threshold == 0.1) {
    height <- 6;
    plot.objects.heights <- c(1.2,
        sapply(peak.clusters, function(cluster) (length(cluster) / length(unlist(peak.clusters)) * 2)),
        2.25);
    } else {
    height <- 4.5;
    plot.objects.heights = c(0.8,
        sapply(peak.clusters, function(cluster) (length(cluster) / length(unlist(peak.clusters)) * 2)),
        1.25);
    }

#using pdf as the font looks strange when saving directly with function 
pdf(width = 6.75, height = height, paste0('results/hypoxia/hypoxia_m6A_peak_summary_FDR_', config$FDR.threshold, 'peak_order_', config$peak.order.method, '_noCNA.pdf'))
create.multipanelplot(
    plots,
    layout.skip = c(FALSE, TRUE, rep(FALSE, length(plots) -1)),
    #filename = paste0('results/hypoxia/hypoxia_m6A_peak_summary_FDR_', config$FDR.threshold, 'peak_order_', config$peak.order.method, '_noCNA.pdf'),
    layout.width = 2,
    layout.height = 4,
    size.units = 'in',
    resolution = 300,
    style = 'Nature',
    plot.objects.heights = plot.objects.heights,
    plot.objects.widths = c(5, 1),
    y.spacing = c(-2, -2.15, -3),
    x.spacing = -2.5,
    ylab.axis.padding = c(-3, 0),
    width = 6.75,
    height = height,
    bottom.padding = -2.5,
    top.padding = 0,
    left.padding = -1.75,
    top.legend.padding = -1,
    right.legend.padding = 0.5,
    right.padding = 0,
    legend = list(
        right = list(
            fun = plot.legend
            )
        ),
    );
dev.off()
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
