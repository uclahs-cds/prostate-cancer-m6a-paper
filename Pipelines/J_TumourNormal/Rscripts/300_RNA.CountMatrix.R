source("000_HEADER.R")



# Preamble ----------------------------------------------------------------
# Compiling a gene by sample matrix for DE

# Loading Data ------------------------------------------------------------

rsem.input = lapply(all_samples, function(sample_id){
  cat(sample_id, "\n")
  tmp = read.delim(file.path("E_RSEM_output", paste0(sample_id, "_Input.genes.results")))[,c("gene_id", "expected_count", "TPM")]
  colnames(tmp)[colnames(tmp) %in% c("expected_count", "TPM")] = paste0(sample_id, ".", colnames(tmp)[colnames(tmp) %in% c("expected_count", "TPM")])
  tmp
})
rsem.input = Reduce(function(x,y) {merge(x,y, by = c("gene_id"), all = T)}, rsem.input)
rownames(rsem.input) = rsem.input$gene_id

# Extract Read Counts
rsem.counts = rsem.input[,grep("expected_count$", colnames(rsem.input))]
colnames(rsem.counts) = gsub(".expected_count", "", colnames(rsem.counts), fixed = T)

# Normalize Data Matrix ---------------------------------------------------

tumour.normal = rsem.counts[,grep("^T|^N", colnames(rsem.counts))]

xenografts = rsem.counts[,grep("^LT", colnames(rsem.counts))]

# Write Output Table ------------------------------------------------------

write.table(
  tumour.normal,
  file = "Tumour.Normal.RSEM.counts.tsv",
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

write.table(
  xenografts,
  file = "Xenografts.RSEM.counts.tsv",
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
