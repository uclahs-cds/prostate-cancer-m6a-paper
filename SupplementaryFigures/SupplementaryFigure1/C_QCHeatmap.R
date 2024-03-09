source("/cluster/home/helenzhu/code/snakemake_M6A_8_PipelineRebuild/Rscripts/000_HEADER.R")



# Preamble ----------------------------------------------------------------
# Making a heatmap of QC metrics

# Loading Data ------------------------------------------------------------

# Read in QC data from file
qc_data = read.delim("Summary_Tables/QC_compiled_for_Hugo.tsv", header = T, check.names = F)

# Positive vs. Negative Metric
signs = read.delim("Summary_Tables/QC_Metrics.Pos.Neg.txt", header = F)
colnames(signs) = c("Metric", "Sign")

# Z-score -----------------------------------------------------------------

# Z Score
plotting_data = apply(qc_data[,2:ncol(qc_data)], 2, function(i) scale(i, center = T, scale = T))
rownames(plotting_data) = qc_data$Sample

# Negative Z scores
neg_z = signs$Metric[which(signs$Sign == "neg")]
for(i in neg_z){
  plotting_data[,i] = -1*plotting_data[,i]
}

# Transpose for visualization, calculating cum(Z) | Z < 0
plotting_data = t(plotting_data)
plotting_data[plotting_data >0 ] <- 0

# Barplot
barplot_data = colSums(plotting_data)
barplot_df = data.frame("Sample" = names(barplot_data), "Sum" = barplot_data, stringsAsFactors = F)
barplot_df = barplot_df[order(-barplot_df$Sum),]
barplot_df$Sample = factor(barplot_df$Sample, levels = barplot_df$Sample)

# Reordering plotting data
plotting_data = plotting_data[,as.character(barplot_df$Sample)]

ylabels = rownames(plotting_data)
ylabels = gsub(".", " ", ylabels, fixed = T)
ylabels = gsub("Percent", "%", ylabels, fixed = T)
ylabels = gsub("percent", "%", ylabels, fixed = T)
ylabels = gsub("%E", "% E", ylabels, fixed = T)

filename = "figures/199_QC.Metrics.Heatmap.pdf"
# pdf(filename, width = 25, height = 13)
# pdf(filename, width = 25, height = 10)
pdf(filename, width = 8, height = 7)

bp = create.barplot(
  Sum ~ Sample,
  barplot_df,
  # Formatting
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.cex = 0, # 0.8,
  yaxis.tck = 0,
  xaxis.tck = 0,
  xlab.label = "",
  ylab.label = "Sum of Z (Z < 0)",
  ylab.cex = 1,
  axes.lwd = 2,
  border.col = "black",
  # Lines
  line.func = function(x) {x = -20},
  line.from = 0,
  line.to = 162,
  line.col = 'darkgrey'
)

hm = create.heatmap(
  plotting_data,
  # Axes labels
  yaxis.lab = ylabels,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.lab = colnames(plotting_data),
  xaxis.cex = 0, # 0.8,
  yaxis.tck = 0,
  xlab.cex = 1,
  xlab.label = "Samples",
  # Clustering
  clustering.method = "none",
  # Colour scheme
  colour.scheme = c("red", "white"),
  colour.centering.value = 0,
  colourkey.labels.at = c(-10:0),
  colourkey.cex = 1,
  at = seq(0, -10, -2),
  same.as.matrix = T,
  # Grid Rows
  grid.row = TRUE,
  row.lines = 26.5-c(6, 12, 18, 22, 24),
  row.colour = "grey",
  row.lwd = 1,
  axes.lwd = 2
)

# Legend
this_legend <- list(
  legend = list(
    colours = c('red', 'white'),
    labels = c('-10','0'),
    title = 'Z score',
    continuous = TRUE,
    height=3,
    pos.x = 0.1
  )
)

legend.grob <- legend.grob(
  legends = this_legend,
  title.just = 'left',
  label.cex = 1,
  title.cex = 1
)

create.multipanelplot(
  plot.objects = list(bp, hm), #  list(metpeak_bp, exomepeak_bp, bp, hm),
  layout.height = 2, # 4
  layout.width = 1,
  plot.objects.heights = c(1, 3), #  c(1, 1, 1.5, 3),
  y.spacing = -1, # c(-2, -2, -1, -2),
  ylab.axis.padding = -18, # 17
  left.padding = 8,
  main = "QC Summary",
  main.cex = 0,
  # legend = list(
  #   right = list(
  #     x = 0.10,
  #     y = 0.50,
  #     fun = legend.grob
  #   )
  # )
)

dev.off()
system(paste0('cp ', filename, ' ~/figures'))

# QQ Plot -----------------------------------------------------------------


filename = "figures/199_QC.Metrics.QQ.Plot.pdf"
pdf(filename, width = 10, height = 10)

df = data.frame(
  "Theoretical" = quantile(runif(10000), probs = seq(0.01,0.99,0.01)),
  "Real" = quantile(barplot_data, probs = seq(0.01,0.99,0.01)),
  stringsAsFactors = F
)

create.scatterplot(
  Real ~ Theoretical,
  df,
  main = "QQ Plot, Uniform"
)

create.histogram(
  barplot_data,
  main = "Histogram",
  xlab.label = "Sum of Negative Z scores",
  ylab.label = "Percent"
)

dev.off()
system(paste0('cp ', filename, ' ~/figures'))

# Checkpoint --------------------------------------------------------------

# Checking what the three worst samples are bad at
plotting_data[,160:162]

which(colnames(plotting_data) == "CPCG0269")
# [1] 143

which(colnames(plotting_data) == "CPCG0587")
#  [1] 160

# Writing out which samples to exclude ------------------------------------

exclude_samples = as.character(barplot_df$Sample[barplot_df$Sum < -20])
save(exclude_samples, file = "Summary_Tables/exclude_samples.rsav")

write.table(
  exclude_samples,
  file = "Summary_Tables/exclude_samples.tsv",
  sep = "\n",
  row.names = F,
  col.names = F,
  quote = F
)
