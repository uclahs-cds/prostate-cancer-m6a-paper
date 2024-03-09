library(exomePeak)

args = commandArgs(trailingOnly=TRUE)

# Initializing
sample_id = args[1]
ip_bam = args[2]
input_bam = args[3]
gtf = args[4]
output_dir = args[5]

# Output Directory
setwd(output_dir)

exomepeak(
  GENE_ANNO_GTF = gtf,
  IP_BAM = ip_bam,
  INPUT_BAM = input_bam,
  EXPERIMENT_NAME = sample_id,
  FRAGMENT_LENGTH=300,
  READ_LENGTH=150)
