source("000_HEADER.R")

library(BoutrosLab.plotting.general)


# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
isamples = args[1]

# QC directory
basewd = "Plink/"

# 133 European samples, filtered by minor allele frequency
# or all 162 samples, unfiltered
# isamples = 133
tag = paste0("CPCGENE_genotype_", isamples)

# SNP-Specific QC ---------------------------------------------------------

# Distribution of minor allele frequency
frq = read.table(
  paste0(basewd, tag, ".frq"),
  header = F,
  sep = "\t",
  skip = 1,
  stringsAsFactors = F)
colnames(frq) = c("CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE:FREQ1", "ALLELE:FREQ2")

filename = "~/figures/CPCGENE_133.MAF.pdf"
pdf(filename, width = 5, height = 5)
create.histogram(
  frq[,'ALLELE:FREQ2'],
  type = "percent",
  breaks = seq(0, 1, 0.05),
  ylab.label = "Percent (%)",
  ylab.cex = 1,
  yaxis.cex = 1,
  xlab.label = "MAF",
  xlab.cex = 1,
  xaxis.cex = 1,
  xlimits = c(0, 1),
  ylimits = c(0, 100),
  main = "133 European Samples",
  main.cex = 1,
  # Adding a line
  abline.v = 0.05,
  abline.col = "red",
  abline.lty = 2
)
dev.off()

# Distribution of Hardy Weinberg equilibrium
hwe = read.table(
  paste0(basewd, tag, ".hwe"),
  header = T,
  sep = "\t",
  stringsAsFactors = F)
hwe.log10 = -log10(hwe[,'P_HWE'])

filename = "~/figures/CPCGENE_133.HWE.pdf"
pdf(filename, width = 5, height = 5)
create.histogram(
  -log10(hwe[,'P_HWE']),
  type = "count",
  breaks = seq(0, 40, 1),
  ylab.label = expression(bold("# of SNPs (10"^"6"*")")),
  ylab.cex = 1,
  yaxis.cex = 1,
  yat = seq(0, 9*10^6, 10^6),
  yaxis.lab = seq(0, 9, 1),
  ylimits = c(0, 9*10^6),
  xlab.label =  expression(bold(paste("-log"["10"]*" HWE P"["value"]))),
  xlab.cex = 1,
  xaxis.cex = 1,
  main = "133 European Samples, MAF > 0.05",
  main.cex = 1,
  # Adding a line
  abline.v = 6,
  abline.col = "red",
  abline.lty = 2
)

create.histogram(
  -log10(hwe[hwe$P_HWE < 10^-6,'P_HWE']),
  type = "count",
  breaks = seq(6, 40, 1),
  ylab.label = expression(bold("# of SNPs (10"^"3"*")")),
  ylab.cex = 1,
  yaxis.cex = 1,
  yat = seq(0, 20000, 5000),
  yaxis.lab = seq(0, 20, 5),
  ylimits = c(0, 20000),
  xlab.label =  expression(bold(paste("-log"["10"]*" HWE P"["value"]))),
  xlab.cex = 1,
  xaxis.cex = 1,
  main = "133 European Samples",
  main.cex = 1,
  # Adding a line
  abline.v = 6,
  abline.col = "red",
  abline.lty = 2
)
dev.off()

# Percent genotyping missing per snp
lmiss = read.table(
  paste0(basewd, tag, ".lmiss"),
  header = T,
  sep = "\t",
  stringsAsFactors = F)

filename = "~/figures/CPCGENE_133.MissingMarkers.pdf"
pdf(filename, width = 5, height = 5)
create.histogram(
  lmiss[,'F_MISS'],
  type = "count",
  breaks = seq(0, 1, 0.05),
  ylab.label = expression(bold("# of SNPs (10"^"6"*")")),
  ylab.cex = 1,
  yaxis.cex = 1,
  yat = seq(0, 10*10^6, 10^6),
  yaxis.lab = seq(0, 10, 1),
  ylimits = c(0, 10*10^6),
  xlab.label =  expression(bold("Missing Marker Proportion")),
  xlab.cex = 1,
  xaxis.cex = 1,
  main = "133 European Samples",
  main.cex = 1,
  # Adding a line
  abline.v = 0.05,
  abline.col = "red",
  abline.lty = 2
)

create.histogram(
  lmiss[lmiss$F_MISS > 0.05,'F_MISS'],
  type = "count",
  breaks = seq(0, 1, 0.05),
  ylab.label = expression(bold("# of SNPs (10"^"4"*")")),
  ylab.cex = 1,
  yaxis.cex = 1,
  yat = seq(0, 6*10^4, 10^4),
  yaxis.lab = seq(0, 6, 1),
  ylimits = c(0, 6*10^4),
  xlab.label =  expression(bold("Missing Marker Proportion")),
  xlab.cex = 1,
  xaxis.cex = 1,
  main = "133 European Samples",
  main.cex = 1,
  # Adding a line
  abline.v = 0.05,
  abline.col = "red",
  abline.lty = 2
)
dev.off()

# Sample-Specific QC ------------------------------------------------------

# Percent genotyping missing in samples
imiss = read.table(
  paste0(basewd, tag, ".imiss"),
  header = T,
  sep = "\t",
  stringsAsFactors = F)

filename = "~/figures/CPCGENE_133.MissingSamples.pdf"
pdf(filename, width = 5, height = 5)

create.histogram(
  imiss[,'F_MISS'],
  type = "count",
  breaks = seq(0, 0.015, 0.001),
  ylab.label = "# Samples",
  ylab.cex = 1,
  yaxis.cex = 1,
  yat = seq(0, 133, 10),
  ylimits = c(0, 135),
  xlab.label = expression(bold("Missing Marker Proportion")),
  xlab.cex = 1,
  xaxis.cex = 1,
  main = "133 European Samples",
  main.cex = 1,
  # Adding a line
  abline.v = 0.05,
  abline.col = "red",
  abline.lty = 2
)

dev.off()


# Writing a table of markers to filter ------------------------------------

lmiss.filter = lmiss[lmiss$F_MISS > 0.05, c("CHR", "POS")]
hwe.filter = hwe[hwe$P_HWE < 10^-6, c("CHR", "POS")]
frq.filter = frq[frq[,'ALLELE:FREQ2'] < 0.05 | frq[,'ALLELE:FREQ1'] < 0.05, c("CHROM", "POS")]
colnames(frq.filter) = c("CHR", "POS")

dim(lmiss.filter)
# [1] 210444      2
dim(hwe.filter)
# [1] 135563      2
dim(frq.filter)
# [1] 2727855       2

filter.variants = rbind(lmiss.filter, hwe.filter, frq.filter)
filter.variants = unique(filter.variants)
filter.variants = filter.variants[complete.cases(filter.variants),]
filter.variants$CHR = factor(filter.variants$CHR, levels = paste0("chr", c(1:22, "X")))
filter.variants = filter.variants[order(filter.variants$CHR, filter.variants$POS),]

dim(filter.variants)
# [1] 2976785       2

# filter.variants = paste0(filter.variants$CHR, ":", filter.variants$POS)
write.table(
  filter.variants,
  file = "Input/Filter.Variants.tsv",
  sep = "\t",
  col.names = F,
  row.names = F,
  quote = F
)

# Extra Metrics -----------------------------------------------------------

# Mean read depth per snp
# ldepth = read.table(
#   paste0(basewd, tag, ".ldepth.mean"),
#   header = T,
#   sep = "\t",
#   stringsAsFactors = F)

# SNP quality
# lqual = read.table(
#   paste0(basewd, tag, ".lqual"),
#   header = T,
#   sep = "\t",
#   stringsAsFactors = F)

# Mean read depth per sample
# idepth = read.table(
#   paste0(basewd, tag, ".idepth"),
#   header = T,
#   sep = "\t",
#   stringsAsFactors = F)
