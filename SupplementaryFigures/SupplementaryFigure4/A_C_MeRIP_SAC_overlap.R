library(GenomicRanges)
library(data.table)
library(ggalluvial)
library(ggplot2)
library(BoutrosLab.plotting.general)

#load SAC-seq peaks
sac.seq <- Sys.glob('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/SAC-Seq/shunliu/outputs/quantification/*.unadjusted.m6A.hits.bed')
names(sac.seq) <- sub('\\..+', '', basename(sac.seq))
sac.seq <- lapply(sac.seq, read.delim, header = FALSE)
sac.seq.GRs <- GRangesList(lapply(sac.seq, makeGRangesFromDataFrame, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V6', keep.extra.columns = TRUE))

#load m6A-seq sample peaks
m6a.seq <- read.delim('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/processed_data/GRCh38/Z_Matrices/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.Mask.tsv')
m6a.seq <- m6a.seq[, grep('CPCG', names(sac.seq), value = T)]
m6a.seq.peaks <- read.delim('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/processed_data/GRCh38/Z_Matrices/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.peaks.tsv', row.names = 1)
#look at overlap of total peak regions analysed

V16 <- makeGRangesFromDataFrame(
    read.delim('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/processed_data/hypoxia_cell_lines/V16A-normoxia/V16A-normoxia_peak.bed6', header = F),
    seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V6'
    )
PC3 <- makeGRangesFromDataFrame(
    read.delim('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/processed_data/hypoxia_cell_lines/PC-3-normoxia/PC-3_normoxia_peak.bed6', header = F),
    seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V6'
    )
m6a.seq.GRs <- GRangesList(c(
    list(V16A = V16, PC3 = PC3),
    apply(m6a.seq, 2, function(sample) makeGRangesFromDataFrame(m6a.seq.peaks[row.names(m6a.seq)[which(as.logical(sample))], ]))
    ))

calc.overlap <- function(sample) {
    intersection.m6A <- length(which(overlapsAny(m6a.seq.GRs[[sample]], sac.seq.GRs[[sample]])))
    intersection.SAC.seq <- length(which(overlapsAny(sac.seq.GRs[[sample]], m6a.seq.GRs[[sample]])))
    data.frame(
        sample = sample,
        m6A.seq = length(m6a.seq.GRs[[sample]]),
        SAC.seq = length(sac.seq.GRs[[sample]]),
        intersection.m6A = intersection.m6A,
        intersection.SAC.seq = intersection.SAC.seq
    )
}

#calculating the overlap between meRIP-seq peaks and SAC-seq sites
full.m6A.GR <- makeGRangesFromDataFrame(m6a.seq.peaks[, 1:6])
peak.data <- do.call(rbind, lapply(names(m6a.seq.GRs), calc.overlap))

#compare all sites
intersection.m6A.all <- length(which(overlapsAny(full.m6A.GR, unique(unlist(sac.seq.GRs)))))
intersection.SAC.seq <- length(which(overlapsAny(unique(unlist(sac.seq.GRs)), full.m6A.GR)))

peak.data <- rbind(data.frame(
    sample = 'Full m6A cohort',
    m6A.seq = length(full.m6A.GR),
    SAC.seq = length(unique(unlist(sac.seq.GRs))),
    intersection.m6A = intersection.m6A.all,
    intersection.SAC.seq = intersection.SAC.seq),
    peak.data)

peak.data$proportion.m6A.seq <- round(peak.data$intersection.m6A / peak.data$m6A.seq, 2)
peak.data$proportion.SAC.seq <- round(peak.data$intersection.SAC.seq / peak.data$SAC.seq, 2)
write.table(
    peak.data,
    'SACseq_m6Aseq_sample_overlap.txt',
    sep = '\t',
    row.names = F,
    col.names = T,
    quote = F
    )

#boxplot of % overlap across samples

to.plot <- data.frame(overlap = peak.data$proportion.SAC.seq[-2] * 100, x = factor(1))
create.boxplot(
    overlap ~ x,
    filename = 'SACseq_m6Aseq_sample_overlap_boxplot.pdf',
    width = 2,
    height = 2.5,
    res  = 300,
    data = to.plot,
    add.stripplot = TRUE,
    xlab.label = '',
    ylab.label = 'Overlapping SAC-seq sites (%)',
    xaxis.lab = '',
    #yaxis.lab = c('',''),
    points.pch = 20,
    alpha = .7,
    main.cex = 1,
    xlab.cex = .7,
    ylab.cex = .7,
    xaxis.cex = .5,
    yaxis.cex = .5,
    left.padding = 2,
    xaxis.tck = 0.2,
    yaxis.tck = 0.3,
    );

#create Sankey diagram detailing overlap between meRIP-seq peaks and SAC-seq sites
SAC.sites <- unique(unlist(sac.seq.GRs[grep('CPCG', names(sac.seq.GRs))]))

tally <- data.frame(
    tech = c(rep('SAC-seq', 2), rep('meRIP-seq', 2)),
    overlap = rep(c(FALSE, TRUE), 2),
    freq = c(table(overlapsAny(SAC.sites, full.m6A.GR)), table(overlapsAny(full.m6A.GR, SAC.sites)))
    )

ggplot(tally,
    aes(y = freq, axis1 = tech, axis2 = overlap)) +
    geom_alluvium(aes(fill = tech), width = 1/12) +
    geom_stratum(width = 1/12, fill = c('#619CFF', '#F8766D', 'black', 'white'), color = "grey") +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Per methodology", "Cross-methodology\n  overlap"), expand = c(.05, .05)) +
    coord_cartesian(xlim = c(1, 1.98)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    geom_text(stat = "flow", aes(label = freq), nudge_x = 0.2, size = 2) +
    ylab(expression(bold(m^6*A~regions))) + 
    theme(
        legend.title=element_blank(), legend.text = element_text(size=6, face = "bold"), legend.key.size = unit(.35, 'cm'),
        axis.title.x = element_blank(), axis.title.y = element_text(colour="black", face="bold", size = 7),
        axis.text.x = element_text(face="bold", colour="black", size = 7),
        axis.text.y = element_text(face="bold", colour="black", size = 6.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "white")
        )
ggsave('full_SACseq_sites_vs_m6A_peaks_alluvium.pdf', height = 3, width = 3, units = 'in')
