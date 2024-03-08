library(exomePeak)

args = commandArgs(trailingOnly=TRUE)

# Initializing
sample_id = args[1]
ip_bam = args[2]
input_bam = args[3]
gtf = args[4]
output_dir = args[5]

# Testing
# sample_id = "CPCG0428"
# ip_bam = "/cluster/home/helenzhu/Cluster_Helen/Snakemake_CPCGENE_m6A/STAR_BAM/CPCG0428_IP_Aligned.sortedByCoord.out.bam"
# input_bam = "/cluster/home/helenzhu/Cluster_Helen/Snakemake_CPCGENE_m6A/STAR_BAM/CPCG0428_Input_Aligned.sortedByCoord.out.bam"
# gtf = "/cluster/home/helenzhu/Cluster_Helen/Snakemake_CPCGENE_m6A/GRCh37_v24lift37/gencode.v24lift37_ek12.gtf"
# output_dir = "/cluster/home/helenzhu/Cluster_Helen/Snakemake_CPCGENE_m6A/MeTPeak_Peaks/"

# Output Directory
setwd(output_dir)

exomepeak(
  GENE_ANNO_GTF = gtf,
  IP_BAM = ip_bam,
  INPUT_BAM = input_bam,
  EXPERIMENT_NAME = sample_id,
  FRAGMENT_LENGTH=300,
  READ_LENGTH=150)
