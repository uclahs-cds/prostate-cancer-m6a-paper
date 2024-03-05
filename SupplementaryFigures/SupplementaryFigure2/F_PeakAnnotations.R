source("/cluster/home/helenzhu/code/snakemake_M6A_9_PeakMotifs/Rscripts/000_HEADER.R")

library(Guitar)
library(valr)

# Preamble ----------------------------------------------------------------
# Makes some plots of the distribution peak annotations across genomic features

# Plots required
# 1. Annotations across gene types
# 2. Distribution across gene structures

# Loading Data ------------------------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

overview.peaks = read.delim(
  "H_PeakCounts/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.bed12",
  header = F
)

# GuitarPlot --------------------------------------------------------------

# GTF
gtf = "GRCh38.p13_Build34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf"
bedfile = "H_PeakCounts/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.bed"

# GuitarPlot
# guitarplot.data = GuitarPlot(
#   txGTF=gtf,
#   pltTxType = "mrna",
#   stBedFiles = bedfile,
#   txGuitarTxdbSaveFile = "GRCh38.p13_Build34.gencode.v34")
# save(guitarplot.data, file = "H_PeakCounts/GuitarPlot.MeTPeak.V9.rsav")
print(load("H_PeakCounts/GuitarPlot.MeTPeak.V9.rsav"))
gp.data = guitarplot.data$data

# Gene annotation boundaries
print(load("GuitarTxdb-GRCh38.p13_Build34.gencode.v34-20220226"))
components = guitarTxdb[['mrna']]$componentWidthAverage_pct
coords = cumsum(components)
text.coords = c(0, coords[1:4]) + diff(c(0, coords))/2

# Loading data from other samples
samples.data = lapply(all_samples, function(i){
  cat(i, "\n")
  print(load(paste0("R_GuitarPlot/MeTPeak_", i, ".rsav")))
  p$data$density
})
samples.data = do.call(cbind, samples.data)
colnames(samples.data) = all_samples

samples.empirical.ci = t(apply(samples.data, 1, function(x)  quantile(x, c(0, 1))))# c(.025, .975)
gp.data = cbind(gp.data, samples.empirical.ci)

# Peak Annotations Plot ---------------------------------------------------

# Calculating proportions of peak annotations
overview.peaks$Gene = gsub("[:].*", "", overview.peaks$V4)
gene_type = merge(overview.peaks, gene_info, by.x = "Gene", by.y = "GeneID")
gene_type = gene_type$Type
gene_type[grepl("pseudogene", gene_type)] <- "pseudogene"
gene_type[grepl("^IG|^Mt|^TEC|^TR", gene_type)] <- "other"

bp.data = data.frame(table(gene_type))
bp.data$PCT = bp.data$Freq/sum(bp.data$Freq)*100
bp.data$Study = "m6A"

# Plotting ----------------------------------------------------------------

# create the simple plot
pg.plt = create.polygonplot(
  formula = NA ~ x,
  data = gp.data,
  max = gp.data[,"0%"],
  min = gp.data[,"100%"],
  # Title
  main = "Peak Density",
  main.just = 'left',
  main.x = 0.03,
  main.y = -0.4,
  main.cex = 2,
  # adjust axes limits
  xlimits = c(0,1),
  ylimits = c(-1,4),
  yat = seq(0, 4, 2),
  ylab.cex = 0,
  xlab.cex = 0,
  yaxis.cex = 2,
  xaxis.cex = 0,
  yaxis.tck = 0,
  xaxis.tck = 0,
  # add fill colour
  col ="dodgerblue",
  border.col = "transparent",
  alpha = 0.1,
  # add middle line
  add.median = TRUE,
  median = gp.data$density,
  median.lty = 1,
  median.lwd = 2,
  # Add rectangles
  add.rectangle = T,
  xleft.rectangle = c(0, coords[1:4]),
  ybottom.rectangle = c(-0.05, -0.1, -0.2, -0.1, -0.05) - 0.2,
  xright.rectangle = coords,
  ytop.rectangle = c(0.05, 0.1, 0.2, 0.1, 0.05) - 0.2,
  col.rectangle = c("black", "darkgrey", "lightgrey", "darkgrey", "black"),
  # alpha.rectangle = 1
  # Add text
  add.text = TRUE,
  text.labels = c("1kbp", "5' UTR", "CDS", "3' UTR", "1kbp"),
  text.x = text.coords,
  text.y = -0.6,
  text.col = 'black',
  text.cex = 1.6,
  text.fontface = 'bold'
)

bp.plt = create.barplot(
  Study ~ PCT,
  data = bp.data,
  # Plot characteristics
  plot.horizontal = T,
  stack = T,
  groups = bp.data$gene_type,
  col = default.colours(12)[9:12],
  border.col = "transparent",
  box.ratio = 20,
  # Axes labels
  main = "Gene Annotations",
  main.cex = 2,
  main.just = 'left',
  main.x = 0.03,
  main.y = -1.2,
  ylab.label = "",
  ylab.cex = 0,
  yaxis.cex = 0,
  yaxis.lab = "",
  xlab.label = "% Peaks",
  xlab.cex = 2,
  xaxis.cex = 2,
  xlimits = c(0, 100),
  xat = seq(0, 100, 20),
  xaxis.tck = 0,
  yaxis.tck = 0
)

covariate.legend <- list(
  legend = list(
    colours = c("dodgerblue", "black"),
    labels = c("Sample Distribution", "Joint Peaks"),
    title = expression(bold(underline('Peak Density'))),
    lwd = 0.5
  ),
  legend = list(
    colours = default.colours(12)[9:12],
    labels = c("lncRNA", "Other", "Protein Coding", "Pseudogene"),
    title = expression(bold(underline('Gene Annotations'))),
    lwd = 0.5
  )
)

side.legend <-legend.grob(
  legends = covariate.legend,
  label.cex = 1.5,
  title.cex = 1.5,
  title.just = 'left',
  title.fontface = 'bold',
  size = 2
)


filename = "~/figures/500_GeneType.pdf"
pdf(filename, width = 10, height = 6)

create.multipanelplot(
  plot.objects = list(pg.plt, bp.plt),
  plot.objects.heights = c(5, 3),
  y.spacing = -5,
  legend = list(
    right = list(
      x = 0.8,
      y = 1,
      fun = side.legend
    )
  ),
  right.legend.padding = 1,
  right.padding = 2
)
dev.off()
