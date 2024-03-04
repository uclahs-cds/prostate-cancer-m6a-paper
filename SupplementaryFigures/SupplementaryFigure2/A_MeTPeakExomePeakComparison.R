source("/cluster/home/helenzhu/code/snakemake_M6A_9_PeakMotifs/Rscripts/000_HEADER.R")

# library(MeRIPtools)
# library(RADAR)


# Preamble ----------------------------------------------------------------
# Making a peaks/sample landscape plot

# Loading Data ------------------------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

metpeak = do.call(rbind, lapply(all_samples, function(this_sample){
  peaks = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/F_MeTPeak_Peaks_SS/", this_sample, "/peak.bed"),
                     header = F, sep = "\t", stringsAsFactors = F)
  peaks$Sample = this_sample
  peaks
}))

exomepeak = do.call(rbind, lapply(all_samples, function(this_sample){
  peaks = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/G_exomePeak_Peaks_SS/", this_sample, "/peak.bed"),
                     header = F, sep = "\t", stringsAsFactors = F)
  peaks$Sample = this_sample
  peaks
}))

# Number of Genes
metpeak_genes = unique(metpeak[,c("V4", "Sample")])
metpeak_genes = data.frame(table(metpeak_genes$Sample), stringsAsFactors = F)

exomepeak_genes = unique(exomepeak[,c("V4", "Sample")])
exomepeak_genes = data.frame(table(exomepeak_genes$Sample), stringsAsFactors = F)

# Number of Peaks
metpeak_peaks = data.frame(table(metpeak$Sample), stringsAsFactors = F)
exomepeak_peaks = data.frame(table(exomepeak$Sample), stringsAsFactors = F)

# Peaks per Gene
peaks_per_gene_metpeak = metpeak
peaks_per_gene_metpeak$peaks_per_gene = paste0(peaks_per_gene_metpeak$Sample, "::", peaks_per_gene_metpeak$V4)
peaks_per_gene_metpeak = data.frame(table(peaks_per_gene_metpeak$peaks_per_gene), stringsAsFactors = F)
peaks_per_gene_metpeak = data.frame(table(peaks_per_gene_metpeak$Freq), stringsAsFactors = F)

peaks_per_gene_exomepeak = exomepeak
peaks_per_gene_exomepeak$peaks_per_gene = paste0(peaks_per_gene_exomepeak$Sample, "::", peaks_per_gene_exomepeak$V4)
peaks_per_gene_exomepeak = data.frame(table(peaks_per_gene_exomepeak$peaks_per_gene), stringsAsFactors = F)
peaks_per_gene_exomepeak = data.frame(table(peaks_per_gene_exomepeak$Freq), stringsAsFactors = F)

# Plotting ----------------------------------------------------------------

# Genes Barplot
exomepeak_genes$Peak = "exomePeak"
metpeak_genes$Peak = "MeTPeak"
genes = rbind(exomepeak_genes, metpeak_genes)
genes$groups = cut(genes$Freq, breaks = seq(0, 8500, 500))
genes$count = 1
dummy = data.frame("Var1" = "", "Freq" = 0, "Peak" = "MeTPeak", "groups" = levels(genes$groups), "count" = 0)
genes = rbind(genes, dummy)
genes_bp = aggregate(count ~ groups + Peak, data = genes, FUN = sum)

mp_genes_bp = create.barplot(
  groups ~ count,
  genes_bp,
  plot.horizontal = T,
  # Groups
  groups = genes_bp$Peak,
  col = default.colours(2),
  # Axes
  xaxis.cex = 1.5,
  # xaxis.rot = 90,
  xlab.label = "Samples",
  yaxis.cex = 0,
  xlab.cex = 1.5,
  ylab.cex = 0,
  yaxis.tck = 0,
  xaxis.tck = 0
)

# Peaks Barplot
metpeak_peaks$Peak = "MeTPeak"
exomepeak_peaks$Peak = "exomePeak"
peaks = rbind(exomepeak_peaks, metpeak_peaks)
peaks$groups = cut(peaks$Freq, breaks = seq(0, 22000, 1000))
peaks$count = 1
dummy = data.frame("Var1" = "", "Freq" = 0, "Peak" = "MeTPeak", "groups" = levels(peaks$groups), "count" = 0)
peaks = rbind(peaks, dummy)
peaks_bp = aggregate(count ~ groups + Peak, data = peaks, FUN = sum)

mp_peaks_bp = create.barplot(
  count ~ groups,
  peaks_bp,
  # Groups
  groups = peaks_bp$Peak,
  col = default.colours(2),
  # Axes
  xaxis.cex = 0,
  xaxis.rot = 90,
  yaxis.cex = 1.5,
  xlab.cex = 0,
  ylab.cex = 1.5,
  ylab.label = "Samples",
  xaxis.tck = 0,
  yaxis.tck = 0,
  yat = seq(0, 30, 10),
  ylimits = c(0, 30)
)

# Scatterplot
metpeak_df = merge(metpeak_peaks[,c("Var1", "Freq")],
                   metpeak_genes[,c("Var1", "Freq")], by = "Var1")
metpeak_df$Peak = "MeTPeak"
exomepeak_df = merge(exomepeak_peaks[,c("Var1", "Freq")],
                     exomepeak_genes[,c("Var1", "Freq")], by = "Var1")
exomepeak_df$Peak = "exomePeak"
peak_df = rbind(metpeak_df, exomepeak_df)

metpeak_scatter = create.scatterplot(
  Freq.y ~ Freq.x,
  data = peak_df,
  # Groups
  groups = peak_df$Peak,
  col = default.colours(2),
  # Labels
  xlab.label = "Peaks",
  ylab.label = "Genes",
  xlab.cex = 1.5,
  ylab.cex = 1.5,
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xlimits = c(0, 22000),
  xat = seq(0, 22000, 10000),
  ylimits = c(0, 8500),
  yat = seq(0, 8500, 2000),
  yaxis.tck = 0,
  xaxis.tck = 0,
  # xaxis.rot = 90,
  key = list(
    text = list(
      lab = c("exomePeak", "MeTPeak"),
      cex = 1.5,
      col = 'black'
    ),
    points = list(
      pch = 19,
      col = default.colours(2),
      cex = 1.5,
      alpha = 1
    ),
    x = 0.05,
    y = 0.95,
    padding.text = 1
  )
)

# Multipanelplot
filename = "figures/110_PeaksLandscape.pdf"
pdf(filename, width = 8.5, height = 8)

create.multipanelplot(
  plot.objects = list(mp_peaks_bp, metpeak_scatter, mp_genes_bp),
  layout.height = 2,
  layout.width = 2,
  layout.skip = c(F, T, F, F),
  plot.objects.heights = c(1, 2),
  plot.objects.widths = c(2, 1),
  y.spacing = -1.5,
  x.spacing = -0.5
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))

# Peaks per gene plot
peaks_per_gene = merge(peaks_per_gene_exomepeak, peaks_per_gene_metpeak, by = "Var1", all = T)
peaks_per_gene[is.na(peaks_per_gene)] <- 0
peaks_per_gene$Var1 = strtoi(as.character(peaks_per_gene$Var1))
peaks_per_gene = peaks_per_gene[order(peaks_per_gene$Var1),]

dummy = colSums(peaks_per_gene[peaks_per_gene$Var1 > 10,])
dummy$Var1 <- ">10"
peaks_per_gene = rbind(peaks_per_gene[peaks_per_gene$Var1 %in% 1:10,], dummy)

peaks_per_gene$Var1 = factor(as.character(peaks_per_gene$Var1), levels = as.character(peaks_per_gene$Var1))
colnames(peaks_per_gene) = c("Var1", "exomePeak", "MeTPeak")
peaks_per_gene$exomePeak = peaks_per_gene$exomePeak/sum(peaks_per_gene$exomePeak)*100
peaks_per_gene$MeTPeak = peaks_per_gene$MeTPeak/sum(peaks_per_gene$MeTPeak)*100
peaks_per_gene = melt(peaks_per_gene)

filename = "figures/110_PeaksPerGene.pdf"
pdf(filename, width = 8.5, height = 4)

create.barplot(
  value ~ Var1,
  peaks_per_gene,
  groups = peaks_per_gene$variable,
  col = default.colours(2),
  # Axes,
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xlab.cex = 1.5,
  ylab.cex = 1.5,
  ylab.label = "% Genes",
  xlab.label = "# Peaks",
  xaxis.tck = 0,
  yaxis.tck = 0
)

dev.off()
system(paste0('cp ', filename, ' ~/figures'))
