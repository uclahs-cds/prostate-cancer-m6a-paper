setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_NormalTumour/")

library(Guitar)

# Preamble ----------------------------------------------------------------
# This script generates the GuitarPlot

# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
meth = args[1]
sample_id = args[2]
gtf = args[3]

if(meth == "MeTPeak"){
  bedfile = paste0("MeTPeak_Peaks/", sample_id, "/peak.bed")
} else if(meth == "exomePeak"){
  bedfile = paste0("exomePeak_Peaks/", sample_id, "/peak.bed")
}

# GTF
# gtf = "GRCh38.p13_Build34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf"

# Filename
filename = paste0("R_GuitarPlot/", meth, "_", sample_id, ".rsav")

# ggtitle
gg_title = paste0(sample_id, " ", meth)

# Making a plot -----------------------------------------------------------

# GuitarPlot
guitarplot = GuitarPlot(txGTF=gtf,
                        pltTxType = "mrna",
                        stBedFiles = bedfile)

# Modifying
p = guitarplot +
  theme_bw() +
  xlab("") +
  ggtitle(gg_title) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(legend.title = element_blank())

save(p, file = filename)
