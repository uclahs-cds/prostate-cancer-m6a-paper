source("000_HEADER.R")

library(vcfR)

# Preamble ----------------------------------------------------------------
# Generating the Input files for MatrixQTL
# 1. SNPs
# 2. Gene Expression - Peak Adjusted Counts
# 3. Covariates
# 4. SNP location
# 5. Gene Location - Peak Location

# all.euro.samples = all.euro.samples[!all.euro.samples %in% exclude_samples]

# SNPs --------------------------------------------------------------------

args = commandArgs(trailingOnly = T)
vcf_file = args[1]
output_tag = args[2]
# vcf_file = "Input/CPCGENE_RiskSNPs.final.vcf.gz"

# Generating SNP matrix
vcf = read.vcfR(vcf_file)
gt.txt = extract.gt(vcf, element = 'GT', as.numeric = FALSE)
gt.vec = structure(c(0, 1, 1, 2), names = c("0/0", "0/1", "1/0", "1/1"))
gt.num = apply(gt.txt, 2, function(x) gt.vec[x])
gt.num = data.frame(gt.num, stringsAsFactors = F)
gt.num$id = rownames(gt.txt)
gt.num = gt.num[,c("id", all.euro.samples)]

# Generating SNP location
snp.loc = getFIX(vcf)
snp.loc = snp.loc[,c("ID", "CHROM", "POS")]
snp.loc[,"ID"] = gt.num$id

snp.loc = data.frame(snp.loc, stringsAsFactors = F)

all(gt.num$id == snp.loc$ID)

# SNP matrix
write.table(
  gt.num,
  file = paste0("Input/", output_tag, ".matrix.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

# SNP location
write.table(
  snp.loc,
  file = paste0("Input/", output_tag, ".location.tsv"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
