source("000_HEADER.R")

library(vcfR)
library(RNOmni)

# Preamble ----------------------------------------------------------------
# Generating the Input files for MatrixQTL
# 1. SNPs
# 2. Gene Expression - Peak Adjusted Counts
# 3. Covariates
# 4. SNP location
# 5. Gene Location - Peak Location

# Dataset Version
args = commandArgs(trailingOnly = T)
meth = args[1]
tag = args[2]

thresh = 0.2

# meth = "MeTPeak"
# tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"

# Expression --------------------------------------------------------------

# Generating Peak Location
peak.loc = read.table(
  paste0("Matrices/", meth, ".", tag, ".peaks.tsv"),
  header = T,
  sep = "\t",
  stringsAsFactors = F)
peak.loc = peak.loc[,c("peak", "chr", "start", "end")]

# Peak Expression
peak.expr = read.table(paste0("Matrices/PeakCounts/", meth, ".", tag, ".PeakCounts.Adjusted.tsv"))
peak.expr$peak = rownames(peak.expr)
peak.expr = peak.expr[peak.loc$peak, c("peak", all.euro.samples)]

all(peak.loc$peak == peak.expr$peak)
# [1] TRUE

# Filtering Data ----------------------------------------------------------

# Filtering for 0's
expr.mat = peak.expr[,all.euro.samples]
nzeros = apply(expr.mat, 1, function(x) sum(x == 0))
keep = which(nzeros < floor(thresh*length(all.euro.samples)))

peak.loc = peak.loc[keep,]
peak.expr = peak.expr[keep, all.euro.samples]

cat(length(keep), " out of ", nrow(expr.mat), " peaks remain after filtering at a threshold of ", thresh, "\n")

# Rank Normalization
peak.expr = t(apply(peak.expr, 1, RNOmni::RankNorm))
peak.expr = data.frame(peak.expr)
peak.expr$peak = rownames(peak.expr)
peak.expr = peak.expr[,c("peak", all.euro.samples)]

all(rownames(peak.expr) == peak.loc$peak)
# [1] TRUE

# Peak Expression
write.table(
  peak.expr,
  file = paste0("Input/", meth, ".", tag, ".Expression.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

# Peak Location
write.table(
  peak.loc,
  file = paste0("Input/", meth, ".", tag, ".Peak.Location.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
