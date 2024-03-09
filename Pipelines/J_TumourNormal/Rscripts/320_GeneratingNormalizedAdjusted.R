source("000_HEADER.R")

library(DESeq2)

# Preamble ----------------------------------------------------------------
# This makes Adjusted Peak Count Matrices

# Useful Functions & Variables --------------------------------------------

saf.colnames = c("Geneid", "Chr", "Start", "End", "Strand", "Length")

clean_names = function(x){
  saf.cols = x[1:6]
  sample.cols = x[7:length(x)]
  sample.cols = gsub(".*BAM", "", sample.cols)
  sample.cols = gsub("IP.*", "", sample.cols)
  sample.cols = gsub("^[.]|[_]$", "", sample.cols)
  names(sample.cols) = NULL
  c(saf.cols, sample.cols)
}

# Loading Data ------------------------------------------------------------

# Tumour Normal
tumour.normal = read.delim("H_PeakCounts/TumourNormal.counts", header = T, comment.char = "#")
colnames(tumour.normal) = clean_names(colnames(tumour.normal))

# Xenograft
xenograft = read.delim("H_PeakCounts/Xenograft.counts", header = T, comment.char = "#")
colnames(xenograft) = clean_names(colnames(xenograft))

# Tumour Normal
tumour.normal.gene = read.delim("Tumour.Normal.RSEM.counts.tsv", header = T)

# Xenograft
xenograft.gene = read.delim("Xenografts.RSEM.counts.tsv", header = T)

# Input Normalization -----------------------------------------------------

normalize.input = function(gene_input_counts){
  input.sizeFactors = DESeq2::estimateSizeFactorsForMatrix(gene_input_counts)
  input.normalized = t(t(gene_input_counts)/input.sizeFactors)
  gene_input_counts
}

tumour.normal.normalized.input = normalize.input(tumour.normal.gene)
xenograft.normalized.input = normalize.input(xenograft.gene)

# Compiling Read Counts ---------------------------------------------------

compile.peak.counts = function(peak.counts){

  # Aggregating Counts
  samps = setdiff(colnames(peak.counts), saf.colnames)
  peak.counts.compiled = aggregate(
    formula = . ~ Geneid,
    data = peak.counts[,c("Geneid", samps)],
    FUN = sum
  )

  # Adjusting matrix
  rownames(peak.counts.compiled) = peak.counts.compiled$Geneid
  peak.counts.compiled = peak.counts.compiled[,samps]

  peak.counts.compiled
}

tumour.normal.peak.counts = compile.peak.counts(tumour.normal)
xenograft.peak.counts = compile.peak.counts(xenograft)

# Calculating Normalized & Adjusted Read Counts ---------------------------

adjust.peak.counts = function(peak.counts.compiled, gene.counts){

  # peak.counts.compiled = tumour.normal.peak.counts
  # gene.counts = tumour.normal.normalized.input

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

  peak.counts.adjusted
}

tumour.normal.adjusted.counts = adjust.peak.counts(
  peak.counts.compiled = tumour.normal.peak.counts,
  gene.counts = tumour.normal.normalized.input
)

xenograft.adjusted.counts = adjust.peak.counts(
  peak.counts.compiled = xenograft.peak.counts,
  gene.counts = xenograft.normalized.input
)

# Writing Output ----------------------------------------------------------

# Adjusted Tumour Normal Counts
write.table(
  tumour.normal.adjusted.counts,
  file = "H_PeakCounts/TumourNormal.Adjusted.Counts.tsv",
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# Adjusted Xenograft Counts
write.table(
  xenograft.adjusted.counts,
  file = "H_PeakCounts/Xenograft.Adjusted.Counts.tsv",
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
