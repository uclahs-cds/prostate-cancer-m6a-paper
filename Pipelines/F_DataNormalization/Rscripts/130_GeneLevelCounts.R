source("000_HEADER.R")



# Preamble ----------------------------------------------------------------
# Compiling gene level counts by summing up peak level counts

# Useful Functions & Variables --------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly = T)
meth = args[1]
tag = args[2]

# meth = "exomePeak"
# tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"

# Loading Peaks Counts
peak.raw.counts = read.table(
  file.path("Z_Matrices", "PeakCounts", paste0(meth, ".", tag, ".PeakCounts.Raw.tsv")),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
peak.raw.counts$Geneid = gsub("[:].*", "", rownames(peak.raw.counts))

# Loading Gene Level Reads
gene.counts = read.table(
  file.path("Z_Matrices", "Normalized.Adjusted_GeneLevel_Counts", "Input.normalized.tsv"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

# Compiling Read Counts ---------------------------------------------------

gene.counts.compiled = aggregate(
  formula = . ~ Geneid,
  data = peak.raw.counts[,c("Geneid", all_samples)],
  FUN = sum
)
rownames(gene.counts.compiled) = gene.counts.compiled$Geneid
gene.counts.compiled = gene.counts.compiled[,all_samples]

# Calculating Normalized & Adjusted Read Counts ---------------------------

# Taking the top 1% of IP read counts in peaks
ave.ip = rowMeans(gene.counts.compiled)
ave.top = rownames(gene.counts.compiled)[order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]]

enrich = as.data.frame(gene.counts.compiled[ave.top,]/gene.counts[ave.top,])
enrich = enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]

size.factors.enrich = DESeq2::estimateSizeFactorsForMatrix(enrich)
gene.counts.normalized = t( t(gene.counts.compiled)/size.factors.enrich )

# Adjusting IP
.noZero = function(x){sapply(x,max,1)}
gene.counts.nozero = t(apply(gene.counts, 1, .noZero))
gene.factor = t( apply(gene.counts.nozero, 1, function(x){x/mean(x)}) )

gene.counts.adjusted = gene.counts.normalized/gene.factor[rownames(gene.counts.normalized),]
gene.counts.adjusted = round(gene.counts.adjusted)

# Writing Output ----------------------------------------------------------

# Raw Gene Counts
write.table(
  gene.counts.compiled,
  file = file.path("Z_Matrices", "PeakSum_GeneCounts", paste0(meth, ".", tag, ".PeakSum_GeneCounts.Raw.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# Normalized Read Counts
write.table(
  gene.counts.normalized,
  file = file.path("Z_Matrices", "PeakSum_GeneCounts", paste0(meth, ".", tag, ".PeakSum_GeneCounts.Normalized.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# Adjusted Read Counts
write.table(
  gene.counts.adjusted,
  file = file.path("Z_Matrices", "PeakSum_GeneCounts", paste0(meth, ".", tag, ".PeakSum_GeneCounts.Adjusted.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
