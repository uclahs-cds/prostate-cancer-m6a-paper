source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")



# Preamble ----------------------------------------------------------------
# Making a summary plot of everything


# Loading Data ------------------------------------------------------------

# Setting variables
meth = "MeTPeak"
tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"
output.tag = 'CPCGENE_filtered.final'
n_components = '0'

# m6A Significant Results
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.rsav")))
# [1] "sig.results" "geno" 

# RNA
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.RNA.rsav")))
# [1] "rna.results"

# Protein
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.protein.rsav")))
# [1] "protein.results"

# Clinical
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.clinical.rsav")))
# [1] "clinical.results"

# BCR
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.bcr.rsav")))
# [1] "bcr.results"

# Analysis ----------------------------------------------------------------

# Fill in NA's
rna.results$bulkRNA.pval[is.na(rna.results$bulkRNA.pval)] <- 1
rna.results$bulkRNA.fdr[is.na(rna.results$bulkRNA.fdr)] <- 1

protein.results$protein.pval[is.na(protein.results$protein.pval)] <-1
protein.results$protein.effectSize[is.na(protein.results$protein.effectSize)] <- 0
protein.results$protein.fdr[is.na(protein.results$protein.fdr)] <- 1

# Extracting signficant results
molecular.results = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("snps", "gene", "GeneName", "Literature"), all = TRUE),
                           list(sig.results, rna.results, protein.results))
molecular.results = molecular.results[molecular.results$bulkRNA.fdr < 0.1 & molecular.results$protein.fdr < 0.1,]

# Clinical results
# clinical = merge(clinical.results, bcr.results, by = "loci")

# All results & Re-correcting for p-values
# all.results = merge(molecular.results, clinical, by.x = "snps", by.y = "loci")
# all.results$ISUP.fdr = p.adjust(all.results$ISUP.p, method = "fdr", n = nrow(all.results))
# all.results$T.fdr = p.adjust(all.results$T.p, method = "fdr", n = nrow(all.results))
# all.results$PSA.fdr = p.adjust(all.results$PSA.p, method = "fdr", n = nrow(all.results))
# all.results$q_value = p.adjust(all.results$p_value, method = "fdr", n = nrow(all.results))
# 
# all.results[all.results$ISUP.fdr < 0.1 | all.results$T.fdr < 0.1 | all.results$PSA.fdr < 0.1 | all.results$q_value < 0.1,]

# Plotting ----------------------------------------------------------------

# Making a dotmap for m6A-RNA-Protein results
filename = "~/figures/360_m6A.QTLs.RNA.protein.pdf"
pdf(filename, width = 4, height = 5)

snp.label = gsub("[:].*", "", molecular.results$snps)
peak.label = gsub(".*[:]", "", molecular.results$gene)
xlabels = paste0(snp.label, " | ", molecular.results$GeneName, " | ", peak.label)
bg.data = -log10(t(molecular.results[,c("bulkRNA.fdr", "protein.fdr")]))
effectSize = t(molecular.results[,c("bulkRNA.effectSize", "protein.effectSize")])

key.func = function(x) { 0.1 + abs(x) ; }
key.sizes = c(2, -1, -0.5, 0, 0.5, 1, 2)

key.tmp = list(
  space = 'right',
  points = list(
    cex =  key.func(key.sizes),
    col = c(rep("dodgerblue2", 3), "transparent", rep("darkorange1", 3)),
    pch = 19
  ),
  text = list(
    lab = as.character(key.sizes),
    cex = 1,
    adj = 1
  ),
  title = expression(beta),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)

create.dotmap(
  #
  x = effectSize,
  bg.data = bg.data,
  colour.scheme = c('white','black'),
  filename = NULL,
  # Scaling
  main.cex = 0,
  main = "",
  xlab.label = "SNP | Gene | Peak",
  xaxis.lab = xlabels,
  xaxis.rot = 90,
  yaxis.lab = c("RNA", "Protein"),
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  ylab.cex = 0,
  xlab.cex = 1,
  xaxis.tck = 0,
  yaxis.tck = 0,
  # Spot size
  spot.size.function = key.func,
  # lwd
  col.lwd = 1, 
  row.lwd = 1,
  lwd = 1,
  # NA
  na.spot.size = 1, 
  na.pch = 4, 
  na.spot.size.colour = 'black',
  # Key
  key = key.tmp,
  key.top = 1,
  # Background colour
  bg.alpha = 1,
  at = c(0, 1, 2, 3),
  # Background colourkey
  colourkey = TRUE,
  colourkey.labels.at = c(0, 1, 2, 3), 
  colourkey.labels = c(
    expression("1"^" "),
    expression("10"^"-1"),
    expression("10"^"-2"),
    expression("10"^"-12")
  ), 
  description = 'Dotmap created by BoutrosLab.plotting.general',
  resolution = 50,
  top.padding = 2,
  bottom.padding = 2
)

dev.off()

