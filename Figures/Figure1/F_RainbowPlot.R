source("/cluster/home/helenzhu/code/snakemake_M6A_9_PeakMotifs/Rscripts/000_HEADER.R")

library(dplyr)
library(StatMeasures)
library(RColorBrewer)


# Preamble ----------------------------------------------------------------
# Creating that plot for genes

# Important Genes ---------------------------------------------------------

prostate_genes = c(
  "KLK3",
  "SPOP",
  "AR",
  "NKX3-1",
  "FOXA1",
  "MYC",
  "TP53",
  "ERG",
  "PTEN",
  "RB1",
  "MED12",
  "CDKNB1",
  "ATM",
  "SCHLAP1",
  "NEAT1",
  "PCAT1",
  "MALAT1")
#   "PCAT19",

m6A_genes = c(
  "WTAP",
  "RMB15",
  "METTL3",
  "METTL14",
  "FTO",
  "ALKBH5",
  "YTHDF1",
  "YTHDF2",
  "YTHDF3",
  "YTHDC1"
)

# Loading Data ------------------------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

peaks.info = read.delim("Z_Matrices/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.peaks.tsv",
                        header = T, stringsAsFactors = F)

peaks.binary = read.delim("Z_Matrices/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.Mask.tsv",
                          header = T, stringsAsFactors = F)

peaks.counts = read.delim("Z_Matrices/PeakCounts/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakCounts.Adjusted.tsv",
                         header = T, stringsAsFactors = F)

# Creating Plotting Data --------------------------------------------------

peaks.binary = peaks.binary[peaks.info$peak,]
peaks.counts = peaks.counts[peaks.info$peak,]

plotting.data = data.frame(
  "peak_id" = peaks.info$peak,
  "n_samples" = rowSums(peaks.binary),
  "median.adjusted.count" = apply(peaks.counts, 1, median),
  "gene_id" = gsub("[:].*", "", peaks.info$peak),
  stringsAsFactors = F
)
plotting.data$gene_name = gene_conversion[plotting.data$gene_id]

# Plotting ----------------------------------------------------------------

# Creating a Barplot
barplot.data = data.frame(table(plotting.data$n_samples))
barplot.data = rbind(data.frame("Var1" = c(-11:0, 149:160), "Freq" = 0),
                     barplot.data)
barplot.data = barplot.data[order(strtoi(as.character(barplot.data$Var1))),]
barplot.data$Var1 = factor(as.character(barplot.data$Var1), levels = as.character(barplot.data$Var1))

peaks_barplot = create.barplot(
  log10(Freq) ~ Var1,
  data = barplot.data,
  # ylab.fontface = 1,
  ylab.label = expression(bold("Peaks")),
  xlab.label = "Samples",
  ylab.cex = 2,
  xlab.cex = 0,
  xaxis.cex = 0,
  yaxis.cex = 1.5,
  ylimits = c(0, 4),
  yaxis.lab = c(expression(bold(paste("10"^"0"))),
                expression(bold(paste("10"^"1"))),
                expression(bold(paste("10"^"2"))),
                expression(bold(paste("10"^"3")))),
  yat = log10(c(1, 10, 100, 1000, 10000)),
  yaxis.tck = 0,
  xaxis.tck = 0,
  border.col = 'black',
  col = "transparent"
)

# Creating a Scatterplot
this_palette = brewer.pal(10, "Spectral")
plotting.data$col = decile(plotting.data$median.adjusted.count)
plotting.data$col = this_palette[plotting.data$col]
plotting.data$log2.median.adjusted.count = log2(plotting.data$median.adjusted.count+1)

# Creating Labels
set.seed(3141)
lab = plotting.data[plotting.data$gene_name %in% c(prostate_genes, m6A_genes),]
lab = lab[order(lab$median.adjusted.count, decreasing = T),]
lab = lab[!duplicated(lab$gene_id),]
lab$degree = sample(c(0, 90, 180, 270), nrow(lab), replace = T)
lab$radius = sample(seq(-0.02, 0.02, 0.005), nrow(lab), replace = T)

peaks_scatterplot = create.scatterplot(
  log2.median.adjusted.count ~ n_samples,
  data = plotting.data,
  # xlimits = c(0, 161),
  ylimits = c(-0.1, 21),
  xlimits = c(-11, 160),
  xat = seq(0, 100, 10)*148/100, # seq(0, 160, 10),
  xaxis.lab = seq(0, 100, 10),
  ylab.cex = 2,
  xlab.cex = 2,
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  ylab.label = expression(bold(paste("log"["2"]*" (Median Adjusted m"^"6"*"A + 1)"))),
  xlab.label = expression(bold("% Samples")),
  yaxis.tck = 0,
  xaxis.tck = 0,
  pch = 1,
  alpha = 0.3,
  col = plotting.data$col,
  xaxis.fontface = "bold",
  yaxis.fontface = "bold",
  xlab.fontface = 1,
  ylab.fontface = 1,
  # Adding Text
  add.text = T,
  text.cex = 1,
  text.fontface = 1,
  text.x = lab$n_samples,
  text.y = lab$log2.median.adjusted.count,
  text.labels = lab$gene_name,
  text.guess.label.position = lab$degree, # ifelse(lab$gene_name %in% c("YTHDF2", "NKX3-1"), 180, 0), # sample(1:360, nrow(lab)), #180,
  text.guess.radius.factor = lab$radius, # sample(seq(0,0.3,0.02), nrow(lab), replace = T), # this_radius_factor,
  text.guess.labels = T,
  # Legend
  key = list(
    text = list(
      lab = as.character(c(1:10)),
      cex = 1.1,
      col = 'black',
      fontface = 'bold'
    ),
    points = list(
      pch = 19,
      col = brewer.pal(10, "Spectral"),
      cex = 1,
      alpha = 1
    ),
    title = expression(bold("Decile")),
    cex = 1,
    x = 0.01,
    y = 0.99,
    padding.text = 3
  )
)

# Legend
covariate.legend <- list(
  legend = list(
    colours = brewer.pal(10, "Spectral"),
    labels = as.character(1:10),
    title = expression(bold(underline("Decile"))),
    pch = 19,
    lwd = 0.5
  )
)

side.legend <-legend.grob(
  legends = rev(covariate.legend),
  label.cex = 1,
  title.cex = 1,
  title.just = 'left',
  title.fontface = 'bold',
  size = 2
)

# Figures -----------------------------------------------------------------

filename = "figures/310_MeTPeakOverviewPlot.pdf"
pdf(filename, width = 15, height = 13.5)
create.multipanelplot(
  plot.objects = list(peaks_barplot, peaks_scatterplot),
  plot.objects.heights = c(1.2, 4),
  plot.objects.widths = 1,
  layout.height = 2,
  layout.width = 1,
  y.spacing = -1,
  top.padding = -2,
  bottom.padding = -2,
  # legend = list(
  #   right = list(
  #     x = 0.7,
  #     y = 1,
  #     fun = side.legend
  #   )
  # )
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))
