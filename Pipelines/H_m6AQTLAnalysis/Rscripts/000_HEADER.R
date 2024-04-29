# Essential Packages
library(BoutrosLab.plotting.general)
library(ggplot2)
library(GenomicRanges)
library(reshape2)

# Working Directory
setwd("<working directory>")

# All samples
# Load all sample ids

# All Euro Samples
# Load European samples

# Gene Name Conversion
print(load("gene_name_conversion.rsav"))

# Batch Info
# Load batch info

# Methods
extract_ensembl_from_bins = function(id_vec){
  unlist(lapply(id_vec, function(i){
    strsplit(i, split = ",")[[1]][1]
  }))
}
