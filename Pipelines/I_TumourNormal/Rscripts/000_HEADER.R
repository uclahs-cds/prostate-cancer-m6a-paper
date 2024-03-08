
# Essential Packages
library(BoutrosLab.plotting.general)
library(ggplot2)
library(GenomicRanges)
library(reshape2)

# Working Directory
setwd("<working directory>")

# Gene Name Conversion
print(load("gene_name_conversion.rsav"))

# All Samples
all_samples = read.table("<path to Tumour.Normal.Config.tsv>",
                         header = T, sep = "\t", stringsAsFactors = F)
all_samples = unique(all_samples$Shortened.SampleID)
