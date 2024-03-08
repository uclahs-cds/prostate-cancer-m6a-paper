source("/cluster/home/helenzhu/code/snakemake_M6A_36_m6AQTLDB/Rscripts/000_HEADER.R")
source("/cluster/home/helenzhu/code/snakemake_M6A_15_TargetIdentification/Rscripts/general_univariate_survival_functions.R")

library(survival)
library(BoutrosLab.statistics.survival)
library(BoutrosLab.plotting.survival)

# Preamble ----------------------------------------------------------------
# This tests for BCR associations for the 36 SNPs

# Loading Data ------------------------------------------------------------

# Testing only the significant snps of allelic imbalance
print(load("Database.m6A.analysis/peak.variants.allelic.imbalance.rsav"))
# [1] "sig.results" "all.results"

# Genotype Data & 36 SNPs
print(load("Database.m6A.analysis/peak.variants.rsav"))
# # [1] "peak.variants" "geno"  

# Survival Data
print(load("SummaryFilesNew/bcr_data.rsav"))
# [1] "bcr_data"
bcr_data$sample = gsub("[-].*", "", rownames(bcr_data))
bcr_data = unique(bcr_data[,c("bcr", "time_to_bcr", "mets", "time_to_mets", "sample")])
bcr_data = bcr_data[bcr_data$sample %in% all.euro.samples,]
rownames(bcr_data) = bcr_data$sample
# bcr_data = bcr_data[complete.cases(bcr_data),]

# Analysis ----------------------------------------------------------------

# Common samples
common.samples = intersect(all.euro.samples, rownames(bcr_data))

# Reformatting
geno.df = geno
rownames(geno.df) = geno.df$id
geno.df = geno.df[sig.results$ID, common.samples]

# Making it binary
geno.df[geno.df == 1] <- 0
geno.df[geno.df == 2] <- 1

# BCR
bcr.results = uni.coxph.all(
  loci.df = geno.df, 
  transpose = TRUE, 
  surv.df = bcr_data, 
  event_col = 'bcr', 
  time_col = 'time_to_bcr', 
  already_binary = TRUE
)

# Metastasis
mets.results = uni.coxph.all(
  loci.df = geno.df, 
  transpose = TRUE, 
  surv.df = bcr_data, 
  event_col = 'mets', 
  time_col = 'time_to_mets', 
  already_binary = TRUE
)

# Updating Significant Results --------------------------------------------

save(bcr.results, mets.results, file = "Database.m6A.analysis/peak.variants.bcr.rsav")


# Plotting ----------------------------------------------------------------

survobj = Surv(bcr_data[all.euro.samples,'time_to_bcr'], bcr_data[all.euro.samples,'bcr']);
patient.groups = factor(as.character(geno.df["rs2240912:chr9:130135457:A:G", all.euro.samples]), levels = c("0", "1"))

filename = "~/figures/150_BCR.m6A.sites.pdf"
pdf(filename, width = 8, height = 10)

create.km.plot(
  survival.object = survobj,
  patient.groups = patient.groups,
  # Axes labels
  xlab.label = expression('Time (years)'),
  xlab.cex = 2,
  ylab.label = expression('BCR-free survival'),
  ylab.cex = 2,
  xaxis.fontface = 1,
  yaxis.fontface = 1,
  main = expression("rs2240912"),
  main.cex = 2,
  # Labels
  risk.labels =  c("AA+AB", "BB"),
  risk.label.fontface = 1,
  risktable.fontsize = 20, 
  key.groups.labels = c("AA+AB", "BB"),
  key.groups.cex = 2,
  key.stats.cex = 2,
  # Padding
  ylab.axis.padding = 1,
  bottom.padding = 2,  
  top.padding = 1,
  right.padding = 1,
  left.padding = 3
)

dev.off()