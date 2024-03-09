source("000_HEADER.R")

# Preamble ----------------------------------------------------------------
# Compiling covariates based on Plink PCA


# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly = T)
input.file = args[1]
output.file = args[2]


# Formatting Input --------------------------------------------------------

covs = read.table(input.file)
rownames(covs) = covs[,1]
covs = covs[,3:ncol(covs)]
covs = t(covs)
rownames(covs) = paste0("Cov", 1:nrow(covs))

write.table(
  covs,
  file = output.file,
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
