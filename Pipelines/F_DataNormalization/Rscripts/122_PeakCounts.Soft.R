source("000_HEADER.R")



# Preamble ----------------------------------------------------------------
# Compiling Peak Counts

# Input files
# 1. peak raw counts (featureCounts)
# 2. gene level counts (RSEM)

# Output files
# 1. raw counts per peak
# 2. peak tpm
# 3. normalized counts per peak
# 4. adjusted counts per peak


# Useful Functions & Variables --------------------------------------------

all_samples = all_samples[!all_samples %in% exclude_samples]

saf.colnames = c("Geneid", "Chr", "Start", "End", "Strand", "Length")

clean_names = function(x){
  saf.cols = x[1:6]
  sample.cols = x[7:length(x)]
  sample.cols = sapply(sample.cols, function(x) strsplit(x, split = "[.]|[_]")[[1]][13])
  names(sample.cols) = NULL
  c(saf.cols, sample.cols)
}

# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly = T)
meth = args[1]
tag = args[2]

# meth = "MeTPeak"
# tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"

# Loading Peaks Counts
peak.counts = read.table(
  file.path("H_PeakCounts", "IP", paste0(meth, ".", tag, ".counts")),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
colnames(peak.counts) = clean_names(colnames(peak.counts))

# Loading Peak Indices
peak.indices = peak.counts[,saf.colnames]

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

# Compiling Read Counts ---------------------------------------------------

peak.counts.compiled = aggregate(
  formula = . ~ Geneid,
  data = peak.counts[,c("Geneid", all_samples)],
  FUN = sum
)

# Calculating TPM ---------------------------------------------------------

# Calculating Peak Length
peak.length = aggregate(
  formula = Length ~ Geneid,
  data = peak.indices,
  FUN = sum
)
peak.length = structure(
  peak.length$Length,
  names = peak.length$Geneid
)

# Peak Counts
rownames(peak.counts.compiled) = peak.counts.compiled$Geneid
peak.counts.compiled = peak.counts.compiled[,all_samples]

# Calculating TPM
peak.counts.tpm = peak.counts.compiled / peak.length[rownames(peak.counts.compiled)]
peak.counts.tpm = t( t(peak.counts.tpm) * 1e6 / colSums(peak.counts.tpm) )

# Calculating Normalized & Adjusted Read Counts ---------------------------

# Taking the top 1% of IP read counts in peaks
ave.ip = rowMeans(peak.counts.compiled)
ave.top = rownames(peak.counts.compiled)[order(ave.ip,decreasing = T)[1:round(0.01*length(ave.ip)[1])]]
ave.top.genes = gsub(":.*", "", ave.top)

enrich = as.data.frame(peak.counts.compiled[ave.top,]/gene.counts[ave.top.genes,])
enrich = enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]

size.factors.enrich = DESeq2::estimateSizeFactorsForMatrix(enrich)
peak.counts.normalized = t( t(peak.counts.compiled)/size.factors.enrich )

# Adjusting IP
.noZero = function(x){sapply(x,max,1)}
gene.counts.nozero = t(apply(gene.counts, 1, .noZero))
gene.factor = t( apply(gene.counts.nozero, 1, function(x){x/mean(x)}) )

peak.genes = gsub(":.*", "", rownames(peak.counts.normalized))
peak.counts.adjusted = peak.counts.normalized/gene.factor[peak.genes,]
peak.counts.adjusted = round(peak.counts.adjusted)


# Masking -----------------------------------------------------------------

mask.mat = function(mat, sig){
  mat = mat[rownames(sig), colnames(sig)]
  cat("Rownames: ", all(rownames(mat) == rownames(sig)), "\n")
  cat("Colnames: ", all(colnames(mat) == colnames(sig)), "\n")
  mat = mat*sig
  mat
}

peak.counts.tpm = mask.mat(mat = peak.counts.tpm, sig = peaks.sig)
peak.counts.normalized = mask.mat(mat = peak.counts.normalized, sig = peaks.sig)
peak.counts.adjusted = mask.mat(mat = peak.counts.adjusted, sig = peaks.sig)


# Writing Output ----------------------------------------------------------

# Peak TPMs
write.table(
  peak.counts.tpm,
  file = file.path("Z_Matrices", "PeakCounts_SoftMask", paste0(meth, ".", tag, ".PeakCounts_SoftMask.TPM.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# Normalized Read Counts
write.table(
  peak.counts.normalized,
  file = file.path("Z_Matrices", "PeakCounts_SoftMask", paste0(meth, ".", tag, ".PeakCounts_SoftMask.Normalized.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# Adjusted Read Counts
write.table(
  peak.counts.adjusted,
  file = file.path("Z_Matrices", "PeakCounts_SoftMask", paste0(meth, ".", tag, ".PeakCounts_SoftMask.Adjusted.tsv")),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
