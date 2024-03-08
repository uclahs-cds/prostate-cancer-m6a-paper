source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")



# Preamble ----------------------------------------------------------------

# Making a summary plot for the risk SNP local QTL analysis


# Loading Data ------------------------------------------------------------

sig.results <- read.delim("results/riskSNP_local_QTL.tsv", header = T)


# Plotting ----------------------------------------------------------------


bg.data = -log10(matrix(sig.results[,"FDR"]))
dot.size = matrix(sig.results[,"beta"])
sig.results$snp_id <- gsub("[:].*", "", sig.results$snps)
sig.results$peak_id <- gsub(".*[:]", "", sig.results$peak)
sig.results$tag <- paste0(sig.results$snp_id, " | ", sig.results$gene_name, " | ", sig.results$peak_id)

filename = "~/figures/514_riskSNP_local_QTL.summary.pdf"
pdf(filename, width = 4, height = 4)

key.func = function(x) { 0.1 + abs(x) ; }
key.sizes = c(-1, -0.5, 0, 0.5, 1)

key.tmp = list(
  space = 'right',
  points = list(
    cex =  key.func(key.sizes),
    col = c(rep("dodgerblue2", 2), "transparent", rep("darkorange1", 2)),
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
  x = dot.size,
  bg.data = bg.data,
  colour.scheme = c('white','black'),
  filename = NULL,
  # Scaling
  main.cex = 0,
  main = "",
  yaxis.lab = sig.results$tag,
  # xaxis.lab = expression(bold("Local Risk SNP m"^"6"*"A QTL")),
  xaxis.lab = "",
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  ylab.label = "SNP | Gene | Peak",
  ylab.cex = 1,
  xlab.cex = 0,
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
    expression("10"^"-3")
  ), 
  description = 'Dotmap created by BoutrosLab.plotting.general',
  resolution = 50
)

dev.off()

