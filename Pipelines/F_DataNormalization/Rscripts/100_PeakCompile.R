source("000_HEADER.R")


# Preamble ----------------------------------------------------------------
# This script collects the data output from the ConsensusPeaks tool and
# writes them into a bed file

# Importing Arguments -----------------------------------------------------

args = commandArgs(trailingOnly = T)
tag = args[1]

cat("Loading Data From Run: ", tag, "\n")

# Loading Data ------------------------------------------------------------

batches = 1:1000

all.methods = c(
  "MeTPeak",
  "exomePeak"
)

peakfiledir = structure(
  c(
    "F_MeTPeak_SF/",
    "G_exomePeak_SF/"
  ),
  names = all.methods
)

load.data = function(meth, tag){

  tmp.file = file.path("ConsensusPeaks", tag, peakfiledir[meth], paste0(tag, ".", meth, ".rsav"))
  if(file.exists(tmp.file)){
    print(load(tmp.file))
  } else {
    peaks = data.frame()
    for(i in batches){
      if(i %% 100 == 0){cat(i, " files loaded\n")}
      filename = file.path("ConsensusPeaks", tag, peakfiledir[meth], paste0(meth, ".", i, ".MergedPeak.tsv"))
      if(!file.exists(filename)){stop(paste0("File ", meth, " ", i, " does not exist!"))}
      tmp = read.delim(filename)
      peaks = rbind(peaks, tmp)
    }
    save(peaks, file = tmp.file)
  }
  return(peaks)
}

# Return a BED12 file
bed12.colnames = c("chr", "start", "end", "name", "score", "strand", "thickStart",
                   "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
retrieve.bed12 = function(peaks){
  # Using the peak label as the gene name
  peaks$name = peaks$peak
  peaks[,bed12.colnames]
}

# Writing Output ----------------------------------------------------------

exomepeak.peaks = load.data(meth = "exomePeak", tag = tag)
exomepeak.bed = retrieve.bed12(exomepeak.peaks)

write.table(
  exomepeak.bed,
  file = file.path("H_PeakCounts", paste0("exomePeak.", tag, ".bed12")),
  col.names = F,
  row.names = F,
  quote = F,
  sep = "\t"
)

write.table(
  exomepeak.peaks,
  file = file.path("Z_Matrices", paste0("exomePeak.", tag, ".peaks.tsv")),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

metpeak.peaks = load.data(meth = "MeTPeak", tag = tag)
metpeak.bed = retrieve.bed12(metpeak.peaks)

write.table(
  metpeak.bed,
  file = file.path("H_PeakCounts", paste0("MeTPeak.", tag, ".bed12")),
  col.names = F,
  row.names = F,
  quote = F,
  sep = "\t"
)

write.table(
  metpeak.peaks,
  file = file.path("Z_Matrices", paste0("MeTPeak.", tag,".peaks.tsv")),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)
