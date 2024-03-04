source("/cluster/home/helenzhu/code/snakemake_M6A_8_PipelineRebuild/Rscripts/000_HEADER.R")

# Preamble ----------------------------------------------------------------
# This calculates the 4 Quality Control Metrics required for MeRIPseq data

# PCR Duplicate Proportion ------------------------------------------------
import_picard_results = function(samplename){
  cat(samplename, "\n")
  tmp = readLines(paste0("D_PicardMarkDup_BAMs/", samplename, ".metrics"))
  idx = grep("^## METRICS", tmp)
  # tmp = tmp[(idx+1):(idx+2)]
  # tmp = do.call(rbind, lapply(tmp, function(i) strsplit(i, split = "\t")[[1]]))
  tmp = tmp[idx+2]
  tmp = strsplit(tmp, split = "\t")[[1]][9]
  data.frame(
    "Sample" = samplename,
    "PDP" = as.numeric(tmp),
    stringsAsFactors = F
  )
}

samples = c(paste0(all_samples, "_IP"), paste0(all_samples, "_Input"))
dup_results = do.call(rbind, lapply(samples, import_picard_results))

paper_stats = dup_results
paper_stats$library = gsub(".*[_]", "", paper_stats$Sample)
summary(paper_stats$PDP[paper_stats$library == "IP"])
summary(paper_stats$PDP[paper_stats$library == "Input"])


# ENCODE Metrics ----------------------------------------------------------

metrics_dir = "D_ENCODE_METRICS/"
import_encode_results = function(samplename){
  cat(samplename, "\n")
  tmp = read.table(paste0(metrics_dir, samplename, ".encode.metrics"), sep = "\t", header = F, stringsAsFactors = F)
  colnames(tmp) = c("TotalReadPairs", "DistinctReadPairs", "OneReadPair", "TwoReadPairs", "NRF.Distinct.Total",  "PBC1.OnePair.Distinct", "PBC2.OnePair.TwoPair")
  tmp$Sample = samplename
  tmp
}
encode_results = do.call(rbind, lapply(samples, import_encode_results))

# Merging Data ------------------------------------------------------------

merged_data = merge(dup_results, encode_results, by = "Sample")
merged_data$Experiment = gsub(".*_", "",  merged_data$Sample)
merged_data$Sample = gsub("_.*", "", merged_data$Sample)

merged_data.ip = merged_data[merged_data$Experiment == "IP",]
colnames(merged_data.ip)[2:ncol(merged_data.ip)] = paste0(colnames(merged_data.ip)[2:ncol(merged_data.ip)], ".IP")

merged_data.input = merged_data[merged_data$Experiment == "Input",]
colnames(merged_data.input)[2:ncol(merged_data.input)] = paste0(colnames(merged_data.input)[2:ncol(merged_data.input)], ".Input")

merged_data = merge(merged_data.ip, merged_data.input, by = "Sample")
merged_data = merged_data[!colnames(merged_data) %in% c("Experiment.Input", "Experiment.IP")]

# PCR Duplicate Proportion
pdp = create.scatterplot(
  PDP.IP ~ PDP.Input,
  data = merged_data,
  main = "PDP",
  main.cex = 1.5,
  main.x = 0.15,
  main.y = -1,
  main.just = "left",
  # xlab.label = "PDP",
  # ylab.label = "PDP",
  xaxis.cex = 1,
  yaxis.cex = 1,
  ylab.cex = 0,
  xlab.cex = 0,
  xat = seq(0, 1, 0.2),
  yat = seq(0, 1, 0.2),
  xlimits = c(0,1),
  ylimits = c(0,1),
  xaxis.tck = 0,
  yaxis.tck = 0
)
# pdp

# Non-Redundant Fraction
nrf = create.scatterplot(
  NRF.Distinct.Total.IP ~ NRF.Distinct.Total.Input,
  data = merged_data,
  main = "NRF",
  main.cex = 1.5,
  main.x = 0.15,
  main.y = -1,
  main.just = "left",
  # xlab.label = "NRF",
  # ylab.label = "NRF",
  xaxis.cex = 1,
  yaxis.cex = 1,
  ylab.cex = 0,
  xlab.cex = 0,
  xat = seq(0, 1, 0.2),
  yat = seq(0, 1, 0.2),
  xlimits = c(0,1),
  ylimits = c(0,1),
  xaxis.tck = 0,
  yaxis.tck = 0
)
# nrf

# PCR Bottlenecking Coefficients #1
pbc1 = create.scatterplot(
  PBC1.OnePair.Distinct.IP ~ PBC1.OnePair.Distinct.Input,
  data = merged_data,
  main = "PBC1",
  main.cex = 1.5,
  main.x = 0.15,
  main.y = -1,
  main.just = "left",
  # xlab.label = "PBC1",
  ylab.label = "PBC1",
  xaxis.cex = 1,
  yaxis.cex = 1,
  ylab.cex = 0,
  xlab.cex = 0,
  xat = seq(0, 1, 0.2),
  yat = seq(0, 1, 0.2),
  xlimits = c(0,1),
  ylimits = c(0,1),
  xaxis.tck = 0,
  yaxis.tck = 0
)

# pbc1

# PCR Bottlenecking Coeffients #2
pbc2 = create.scatterplot(
  PBC2.OnePair.TwoPair.IP ~ PBC2.OnePair.TwoPair.Input,
  data = merged_data,
  main = "PBC2",
  main.cex = 1.5,
  main.x = 0.1,
  main.y = -1,
  main.just = "left",
  # xlab.label = "PBC2",
  # ylab.label = "PBC2",
  xaxis.cex = 1,
  yaxis.cex = 1,
  ylab.cex = 0,
  xlab.cex = 0,
  xat = seq(0, 10, 2),
  yat = seq(0, 10, 2),
  xlimits = c(0,10),
  ylimits = c(0,10),
  xaxis.tck = 0,
  yaxis.tck = 0
)
# pbc2

filename = paste0("figures/170_ENCODE_QC.pdf")
pdf(filename, width = 12, height = 3.8)

create.multipanelplot(
  plot.objects=list(pdp, nrf, pbc1, pbc2),
  # plot.objects.heights = rep(1, 4),
  plot.objects.widths = c(1, 1, 1, 1),
  layout.width = 4,
  layout.height = 1,
  x.spacing = 0.5,
  y.spacing = -0.5,
  ylab.cex = 1.5,
  ylab.label = "IP Library",
  xlab.cex = 1.5,
  xlab.label = "Input Library",
  ylab.axis.padding = 0, # c(-4, -4, -4, -4),
  xlab.axis.padding = 0, # c(-4, -4, -4, -4),
  top.padding = 0,
  left.padding = 1,
  bottom.padding = 1,
  right.padding = 1,
  left.legend.padding = 0,
  right.legend.padding = 0, 
  bottom.legend.padding = 0, 
  top.legend.padding = 1,
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))

# Write a Table for Supplementary -----------------------------------------

merged_data = merged_data[,c("Sample",
                             "PDP.IP", "PDP.Input",
                             "NRF.Distinct.Total.IP", "NRF.Distinct.Total.Input",
                             "PBC1.OnePair.Distinct.IP",  "PBC1.OnePair.Distinct.Input",
                             "PBC2.OnePair.TwoPair.IP",  "PBC2.OnePair.TwoPair.Input" )]
write.table(
  merged_data,
  file = "SummaryFilesNew/ENCODE.QC.tsv",
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)
