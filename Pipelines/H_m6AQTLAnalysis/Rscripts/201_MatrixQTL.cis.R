source("000_HEADER.R")

library("MatrixEQTL")

# Preamble ----------------------------------------------------------------
# Trying to run MatrixQTL for a basic m6a-qtl analysis

# Data --------------------------------------------------------------------

# Version of analysis
args = commandArgs(trailingOnly = TRUE)
meth = args[1]
tag = args[2]
output_tag = args[3]
n_components = args[4]

# Basic Analysis
base.dir = "Input/"
SNP_file_name = paste(base.dir, output_tag, ".matrix.tsv", sep="")
expression_file_name = paste(base.dir, meth, ".", tag, ".Expression.tsv", sep="")
covariates_file_name = paste0("covariates/", meth, ".", tag, "_", n_components, "_peer.tsv")

# Cis & trans QTLs
snps_location_file_name = paste(base.dir, output_tag, ".location.tsv", sep="")
gene_location_file_name = paste(base.dir, meth, ".", tag, ".Peak.Location.tsv", sep="")

# SNPs file
snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

# Gene expression file
gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

# Covariates file
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

# Output Files ------------------------------------------------------------

output_file_name = paste0("results/",  meth, ".", tag, ".", output_tag, ".", n_components, ".MatrixQTL.tsv")
output_file_name_cis = paste0("results/",  meth, ".", tag, ".", output_tag, ".", n_components, ".MatrixQTL.cis.all.tsv")
output_file_name_tra = paste0("results/",  meth, ".", tag, ".", output_tag, ".", n_components, ".MatrixQTL.trans.all.tsv")

# Model -------------------------------------------------------------------

# Specify model
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS

# p-value threshold
pvOutputThreshold = 1e-2
pvOutputThreshold_cis = 1
pvOutputThreshold_tra = if(output_tag == 'CPCGENE_RiskSNPs.final') 1 else 0 # 1e-2

# Distance between SNP and gene
cisDist = if(output_tag == 'CPCGENE_RiskSNPs.final') 10e9 else 1e4

# Covariance
errorCovariance = numeric()

# Analysis ----------------------------------------------------------------

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

cat("finished!")

save(me, file = paste0("results/",  meth, ".", tag, ".", output_tag, ".", n_components, ".all.rsav"))
