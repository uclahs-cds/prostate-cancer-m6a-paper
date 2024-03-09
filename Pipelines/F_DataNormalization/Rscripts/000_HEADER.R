# Essential Packages
library(BoutrosLab.plotting.general)
library(ggplot2)
library(GenomicRanges)
library(reshape2)

# Working Directory
setwd("")

# All samples
# Load all samples

# Methods
extract_ensembl_from_bins = function(id_vec){
  unlist(lapply(id_vec, function(i){
    strsplit(i, split = ",")[[1]][1]
  }))
}
