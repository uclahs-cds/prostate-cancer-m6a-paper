
# Essential Packages
library(BoutrosLab.plotting.general)
library(ggplot2)
library(GenomicRanges)
library(reshape2)

# Working Directory
setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38")


# All samples
print(load("/cluster/home/helenzhu/Cluster_Helen/Snakemake_CPCGENE_m6A/Summary_Tables/all_samples.rsav"))

# Gene Name Conversion
print(load("gene_name_conversion.rsav"))

# Exclude Samples
exclude_samples = c("CPCG0587", "CPCG0269")

# Batch Info
batch_info = read.table("Summary_Tables/CPCGENE_BatchInfo.tsv",
                        header = T, sep = "\t", stringsAsFactors = F)
batch_info = batch_info[order(batch_info$Sample),]


# Methods
extract_ensembl_from_bins = function(id_vec){
  unlist(lapply(id_vec, function(i){
    strsplit(i, split = ",")[[1]][1]
  }))
}
