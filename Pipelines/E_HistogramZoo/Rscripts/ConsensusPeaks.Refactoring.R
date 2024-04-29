
library(ConsensusPeaks)

# Preamble ----------------------------------------------------------------
# This runs ConsensusPeaks for peak merging exomePeak and MeTPeak output

# Arguments ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
i = strtoi(args[1])
meth = args[2]
gtf = args[3]
tag = args[4]

# Loading Data ------------------------------------------------------------

# exomePeak & MeTPeak Genes
if(meth == "exomePeak"){
  outputdir = paste0("ConsensusPeaks/", tag, "/G_exomePeak_SF/")
  peakfiledir = "G_exomePeak_Peaks_SS/"
} else if (meth == "MeTPeak"){
  outputdir = paste0("ConsensusPeaks/", tag, "/F_MeTPeak_SF/")
  peakfiledir = "F_MeTPeak_Peaks_SS/"
}

# genes
genes = read.table(paste0(meth, ".genes.ss.tsv"), header = T, sep = "\t", stringsAsFactors = F)
genes = genes$gene[genes$batch == i]

# output file tag
tag = paste0(meth, ".", i)

# Useful Genes
# print(load("useful.genes.rsav"))

# Loading Peaks -----------------------------------------------------------

all_samples = setdiff(all_samples, exclude_samples)
filenames = paste0(peakfiledir, all_samples, "/peak.bed")

list.hist = transcript.bed.to.hist(
  filenames = filenames,
  n_fields = 12,
  gtf.file = gtf,
  gene.or.transcript = "gene",
  histogram.bin.size = 1,
  regions.of.interest = genes,
  comment = "#"
)

# Running Tool ------------------------------------------------------------

ftc.res = bulk.segment.fit(
  coverage.model.obj = list.hist,
  histogram.count.threshold = 0,
  eps = 0.005,
  seed = 123,
  truncated.models = TRUE,
  uniform.peak.threshold = 0.8,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  min.gap.size = 0,
  min.peak.size = 100,
  max.uniform = F,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
)

# Formatting Results ------------------------------------------------------

res = summarize.results(
  coverage.model.obj = ftc.res,
  output.format = "BED12"
)

# Writing Results ---------------------------------------------------------

filename = file.path(outputdir, paste0(meth, ".", i, ".MergedPeak.tsv"))
write.table(
  res,
  file = filename,
  quote = F,
  sep = "\t",
  col.names = T,
  row.names = F
)
