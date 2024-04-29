source("000_HEADER.R")

library(valr)
library(GenomicRanges)

# Preamble ----------------------------------------------------------------
# This creates a SAF file by taking the union of the peaks using ConsensusPeaks
# Identifying joint peaks using the union method for each dataset

# Functions ---------------------------------------------------------------

base0.to.base1 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)+1
  return(gr)
}

base1.to.base0 = function(gr){
  GenomicRanges::start(gr) = GenomicRanges::start(gr)-1
  return(gr)
}

bed6tobed12 = function(
  gr
){
  bed12.gr = range(gr)
  # blockCount
  S4Vectors::mcols(bed12.gr)$blockCount = length(gr)
  # blockSizes
  blockSizes = IRanges::width(gr) - 1
  S4Vectors::mcols(bed12.gr)$blockSizes = paste0(blockSizes, ",", collapse = "")
  # blockStarts
  blockStarts = GenomicRanges::start(gr) - GenomicRanges::start(bed12.gr)
  S4Vectors::mcols(bed12.gr)$blockStarts = paste0(blockStarts, ",", collapse = "")
  # Return GRanges object
  bed12.gr
}

# Loading Data ------------------------------------------------------------

# tumour.normal = grep("^N|^T", all_samples, value = T)
# peak.files = file.path("MeTPeak_Peaks", tumour.normal, "peak.bed")

xenograft = grep("^LT", all_samples, value = T)
peak.files = file.path("MeTPeak_Peaks", xenograft, "peak.bed")

# Loading Files as BED12
# peaks = read_bed(filename = peak.files, n_fields = 12, comment = "#")
peaks = do.call(rbind, lapply(peak.files, function(i) read.delim(i, sep = "\t", header = F, comment.char = "#", stringsAsFactors = F)))
colnames(peaks) = names(valr:::bed12_coltypes)
peaks.gr = makeGRangesFromDataFrame( peaks, keep.extra.columns = T)
peaks.gr = base0.to.base1(peaks.gr)

# Creating BED6 Files
peaks.bed6 = bed12_to_exons( peaks )
peaks.bed6.gr = makeGRangesFromDataFrame( peaks.bed6, keep.extra.columns = T )
peaks.bed6.gr = base0.to.base1(peaks.bed6.gr)

# Creating Merged Peaks ---------------------------------------------------

# Split peaks by gene and reduce them into the union of these peaks
peaks.gr = split(peaks.gr, peaks.gr$name)
peaks.gr = reduce(peaks.gr)
peaks.gr = lapply(seq_along(peaks.gr), function(i) {
  tmp.gr = peaks.gr[[i]]
  mcols(tmp.gr)$gene = names(peaks.gr)[i]
  mcols(tmp.gr)$peak = paste0(names(peaks.gr)[i], "::", 1:length(tmp.gr))
  tmp.gr
})
peaks.gr = do.call(c, peaks.gr)
peaks.gr = split(peaks.gr, peaks.gr$gene)

# Split peaks by gene in BED6 format and reduce them into the union of these peaks
peaks.bed6.gr = split(peaks.bed6.gr, peaks.bed6.gr$name)
peaks.bed6.gr = reduce(peaks.bed6.gr)

# Matching split peaks
gene.names = sort(unique(names(peaks.gr)))
annotate.bed6 = lapply(gene.names, function(gene.name){
  cat(gene.name, "\n")
  bed12 = peaks.gr[[gene.name]]
  bed6 = peaks.bed6.gr[[gene.name]]
  ovl = findOverlaps(bed6, bed12)
  # Safety check
  if( !all(queryHits(ovl) == seq_along(bed6)) ){
    stop("Something went wrong with ", gene.name)
  }
  mcols(bed6)$gene = gene.name
  mcols(bed6)$peak = mcols(bed12)$peak[subjectHits(ovl)]
  bed6
})
annotate.bed6 = do.call(c, annotate.bed6)

# Creating a SAF File -----------------------------------------------------

# No Scientific Notation in Output
options(scipen = 999)

bed6.data.frame = data.frame(annotate.bed6)

# SAF format (This is a base 1 format)
saf = bed6.data.frame[,c("peak", "seqnames", "start", "end", "strand")]
colnames(saf) = c("GeneID", "Chr", "Start", "End", "Strand")

# Output
# filename = "H_PeakCounts/TumourNormal.peaks.saf"
filename = "H_PeakCounts/Xenograft.peaks.saf"
write.table(
  saf,
  file = filename,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)


# Creating a BED12 file ---------------------------------------------------

# Changing back to base 0, splitting by peak, created a bed 12 annotation per peak
annotate.bed6 = base1.to.base0(annotate.bed6)
annotate.bed6 = split(annotate.bed6, annotate.bed6$peak)
annotate.bed12 = lapply(seq_along(annotate.bed6), function(i) {
  tmp = bed6tobed12(annotate.bed6[[i]])
  tmp$peak = names(annotate.bed6)[i]
  tmp
})
annotate.bed12 = do.call(c, annotate.bed12)
annotate.bed12 = data.frame(annotate.bed12)

# Adding extra columns
annotate.bed12$score = "."
annotate.bed12$thickStart = annotate.bed12$start
annotate.bed12$thickEnd = annotate.bed12$end
annotate.bed12$itemRgb = 0

# Formatting
annotate.bed12 = annotate.bed12[,c("seqnames", "start", "end", "peak", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")]

# Writing an output file
# filename = "H_PeakCounts/TumourNormal.peaks.bed"
filename = "H_PeakCounts/Xenograft.peaks.bed"
write.table(
  annotate.bed12,
  file = filename,
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F
)
