source("000_HEADER.R")


# Preamble ----------------------------------------------------------------
# Compiling gene level counts by summing up peak level counts
# Trying to create soft masked gene level counts using Rupert's strategy
#

# Useful Functions & Variables --------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly = T)
meth = args[1]
tag = args[2]

# meth = "MeTPeak"
# tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"

# Loading Peaks Counts
peak.counts.compiled = read.table(
  file.path("Z_Matrices", "PeakCounts", paste0(meth, ".", tag, ".PeakCounts.Raw.tsv")),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
# peak.raw.counts$Geneid = gsub("[:].*", "", rownames(peak.raw.counts))

# Loading Gene Level Reads
gene.counts = read.table(
  file.path("Z_Matrices", "Normalized.Adjusted_GeneLevel_Counts", "Input.normalized.tsv"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

# Loading Significance Threshold
peaks.sig = read.table(
  file.path("Z_Matrices", paste0(meth, ".", tag, ".peaks.tsv")),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
rownames(peaks.sig) = peaks.sig$peak
peaks.sig = peaks.sig[,all_samples]

# Creating a Mask ---------------------------------------------------------

peaks.sig[peaks.sig < 0.05] <- 10
peaks.sig[peaks.sig <= 1] <- 0
peaks.sig[peaks.sig == 10] <- 1

# Estimating Normalization Factors ----------------------------------------

# Taking the top 1% of IP read counts in peaks
ave.ip = rowMeans(peak.counts.compiled)
ave.top = rownames(peak.counts.compiled)[order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]]
ave.top.genes = gsub(":.*", "", ave.top)

enrich = as.data.frame(peak.counts.compiled[ave.top,]/gene.counts[ave.top.genes,])
enrich = enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]

size.factors.enrich = DESeq2::estimateSizeFactorsForMatrix(enrich)

# Masking -----------------------------------------------------------------

peak.counts.compiled = peak.counts.compiled[rownames(peaks.sig),]
all(rownames(peak.counts.compiled) == rownames(peaks.sig))
all(colnames(peak.counts.compiled) == colnames(peaks.sig))
peak.counts.compiled = peak.counts.compiled*peaks.sig

# Compiling Read Counts ---------------------------------------------------

peak.counts.compiled$Geneid = gsub("[:].*", "", rownames(peak.counts.compiled))

gene.counts.compiled = aggregate(
  formula = . ~ Geneid,
  data = peak.counts.compiled[,c("Geneid", all_samples)],
  FUN = sum
)
rownames(gene.counts.compiled) = gene.counts.compiled$Geneid
gene.counts.compiled = gene.counts.compiled[,all_samples]

# Calculating Normalized & Adjusted Read Counts ---------------------------

# Normalizing gene counts using pre-estimated size factors
gene.counts.normalized = t( t(gene.counts.compiled)/size.factors.enrich )

# Adjusting IP
.noZero = function(x){sapply(x,max,1)}
gene.counts.nozero = t(apply(gene.counts, 1, .noZero))
gene.factor = t( apply(gene.counts.nozero, 1, function(x){x/mean(x)}) )

gene.counts.adjusted = gene.counts.normalized/gene.factor[rownames(gene.counts.normalized),]
gene.counts.adjusted = round(gene.counts.adjusted)

# Writing Output ----------------------------------------------------------

# Normalized Read Counts
write.table(
  gene.counts.normalized,
  file = file.path("Z_Matrices", "PeakSum_GeneCounts_Soft", paste0(meth, ".", tag, ".PeakSum_GeneCounts_Soft.Normalized.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# Adjusted Read Counts
write.table(
  gene.counts.adjusted,
  file = file.path("Z_Matrices", "PeakSum_GeneCounts_Soft", paste0(meth, ".", tag, ".PeakSum_GeneCounts_Soft.Adjusted.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
