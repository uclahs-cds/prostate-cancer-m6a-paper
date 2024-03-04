source("/cluster/home/helenzhu/code/snakemake_M6A_8_PipelineRebuild/Rscripts/000_HEADER.R")




# Preamble ----------------------------------------------------------------
# Recreating Julie's correlation histograms



# Loading Data ------------------------------------------------------------

input_rnaseq_data <- "Z_Counts_compiled/2021-05-05_RSEM_Gene_TPM_for_Julie.tsv"
input_rnaseq <- read.delim(input_rnaseq_data, header = T)
colnames(input_rnaseq) <- gsub("[_].*", "", colnames(input_rnaseq))

bulk_rnaseq_data <- "Summary_Tables/Chen.2019.hg38/2020-10-15_RSEM_gene_TPM.txt"
bulk_rnaseq <- read.delim(bulk_rnaseq_data, header = T)
colnames(bulk_rnaseq) <- gsub("[.].*", "", colnames(bulk_rnaseq))

dim(input_rnaseq)
# [1] 67071   162
dim(bulk_rnaseq)
# [1] 67071   144
length(intersect(rownames(input_rnaseq), rownames(bulk_rnaseq)))
# [1] 67071

# Common samples
common_samples <- intersect(colnames(input_rnaseq), colnames(bulk_rnaseq))
input_rnaseq <- input_rnaseq[,common_samples]
bulk_rnaseq <- bulk_rnaseq[,common_samples]

# Genes
genes <- rownames(bulk_rnaseq)

# Analysis ----------------------------------------------------------------

# Calculating gene-wise correlations
genewise_correlations <- unlist(lapply(genes, function(i){
  cor(
    as.numeric(t(bulk_rnaseq[i,])), 
    as.numeric(t(input_rnaseq[i,])), 
    method = "spearman"
  )
}))

# Calculating sample-wise correlations
# 105 samples
samplewise_correlations <- unlist(lapply(common_samples, function(i){
  cor(
    as.numeric(bulk_rnaseq[,i]), 
    as.numeric(input_rnaseq[,i]), 
    method = "spearman"
  )
}))

# Creating histograms
filename <- "~/figures/181_input_bulk_correlations.pdf"
pdf(filename, width = 3, height = 3)

create.histogram(
  samplewise_correlations,
  # Formatting
  ylab.label = "Samples",
  xlab.label = expression(bold("Spearman's "*rho)),
  xlab.cex = 1,
  ylab.cex = 1,
  xaxis.cex = 1,
  yaxis.cex = 1,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xlimits = c(0.7, 0.9),
  xat = seq(0.7, 0.9, 0.05)
)


create.histogram(
  genewise_correlations[!is.na(genewise_correlations)],
  # Formatting
  ylab.label = "Genes",
  xlab.label = expression(bold("Spearman's "*rho)),
  xlab.cex = 1,
  ylab.cex = 1,
  xaxis.cex = 1,
  yaxis.cex = 1,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xlimits = c(-0.5, 1),
  xat = seq(-0.5, 1, 0.5)
)

dev.off()

# Saving data so we won't have to rerun
save(genewise_correlations, samplewise_correlations, file = "bulk_input_correlations.rsav")
