source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")

library("MatrixEQTL")


# Preamble ----------------------------------------------------------------
# Making a few preliminary plots for m6a-QTL analysis

# Loading Data ------------------------------------------------------------

# ['CPCGENE_filtered.final', 'CPCGENE_RiskSNPs.final', 'CPCGENE_50_5_0.5', 'CPCGENE_50_5_0.8', 'CPCGENE_50_5_0.2']

# Setting variables
meth = "MeTPeak"
tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"
output.tag = 'CPCGENE_filtered.final'
n_components = '0'

# Saved Result File
print(load(paste0("results/", meth,  ".", tag, ".", output.tag, ".", n_components, ".all.rsav")))
# [1] "me"
cis.qtl = me$cis$eqtls

# Loading SNP Positions
snp.loc.file = paste0("Input/", output.tag, ".location.tsv")
snp.loc = read.table(
  file = snp.loc.file,
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

# Include Labels
print(load( paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".labels.rsav") ))
# [1] "labels"

# Organizing Data ---------------------------------------------------------

# Adding Position
plotting.data = merge(snp.loc, cis.qtl, by.x = "ID", by.y = "snps", all.y = T)

# Ordering data
chr = paste0("chr", c(1:22, "X", "Y"))
plotting.data$CHROM = factor(plotting.data$CHROM, levels = chr)
plotting.data = plotting.data[order(plotting.data$CHROM, plotting.data$POS),]
plotting.data$ind = 1:nrow(plotting.data)

# Label formatting
plotting.data.tag = paste0(plotting.data$ID, ":", plotting.data$gene)
labels.tag = paste0(labels$snp, ":", labels$gene)
labels$ind = plotting.data$ind[match(labels.tag, plotting.data.tag)]

# Plotting Things
# colours
chr.col = gsub("chr", "", plotting.data$CHROM)
chr.colours = force.colour.scheme(chr.col, scheme = "chromosome")
chr.colours = c(chr.colours, rep('black', nrow(labels)))
# pch
label.pch.scheme = structure(c(15, 15, 18, 20), names = c("m6A.Site.TagSNP" , "m6A.Site", "Literature.TagSNP", "Top.Result"))
pch.all = c(rep(20, nrow(plotting.data)), label.pch.scheme[labels$Label])
# alpha
alpha.all = c(rep(0.1, nrow(plotting.data)), rep(1, nrow(labels)))

# make chr covariate and chr labels
chr.n.genes       <- vector();
chr.tck           <- vector();
chr.pos.genes     <- vector();
chr.break         <- vector();
chr.break[1]      <- 0;

# get a list of chromosomes to loop
chr = paste0("chr", c(1:22, "X", "Y"));

# loop over each chromosome
for ( i in 1:length(chr) ) {

  # get the number of genes that belong to one chromosome
  n <- sum(plotting.data$CHROM == chr[i]);

  # calculate where the labels go
  chr.n.genes[i]   <- n;
  chr.break[i+1]   <- n + chr.break[i];
  chr.pos.genes[i] <- floor(chr.n.genes[i]/2);
  chr.tck[i]       <- chr.pos.genes[i] + which(plotting.data$CHROM == chr[i])[1];
}

# Updating plotting data
plotting.data.simple = rbind(plotting.data[,c("pvalue", "ind")], labels[,c("pvalue", "ind")])

# Updating labels
prostate.cancer.genes = c("PMS2", "KLK2", "MAX", "ETV6", "KMT2C",  "FOXA1" , "SMAD3", "PTCH1",  "PIK3R1", "CBR3-AS1")
m6a.site.genes = c("PDCD11", "CACNA1H", "IL7R", "SPINT1", "LZTR1", "NMT-MINDY4")
other = c("AHNAK2", "RTL6", "TBC1D16", "HLA-A", "ANKEF1")
labels.txt = labels[labels$Label == "Top.Result" | labels$GeneName %in% c(prostate.cancer.genes, m6a.site.genes, other),]

filename = paste0("figures/", meth, ".", tag, ".", output.tag, ".", n_components, ".Manhattan.png")
create.manhattanplot(
  #
  formula = -log10(pvalue) ~ ind,
  data = plotting.data.simple,
  # Saving file
  filename = filename,
  height = 8,
  width = 12,
  resolution = 500,
  # Labels
  xlab.label = expression(bold('Chromosomes')),
  ylab.label = expression(bold("-log"["10"]*" P"["value"])),
  xat = chr.tck,
  xaxis.lab = c(1:22, 'X', 'Y'),
  xaxis.tck = 0,
  xaxis.cex = 0.7,
  yaxis.cex = 0.7,
  yat = seq(0,30,3),
  ylab.cex = 1,
  xlab.cex = 1,
  yaxis.lab = c(
    expression(bold("1")),
    expression(bold("10"^"-3")),
    expression(bold("10"^"-6")),
    expression(bold("10"^"-9")),
    expression(bold("10"^"-12")),
    expression(bold("10"^"-15")),
    expression(bold("10"^"-18")),
    expression(bold("10"^"-21")),
    expression(bold("10"^"-24")),
    expression(bold("10"^"-27")),
    expression(bold("10"^"-30"))
  ),
  xaxis.fontface = "bold",
  yaxis.fontface = "bold",
  alpha = alpha.all,
  col = chr.colours,
  pch = pch.all,
  cex = 0.75,
  # draw horizontal line
  # abline.h = -log10(0.05),
  # abline.lty = 2,
  # abline.lwd = 1,
  # abline.col = 'black',
  add.text = T,
  text.labels = labels.txt$GeneName,
  text.x = labels.txt$ind,
  text.y = -log10(labels.txt$pvalue)+0.5,
  text.col = 'black',
  text.cex = 0.75,
  text.fontface = 1,
  # Adding a legend
  key = list(
    text = list(
      lab = c(expression(underline("Annotation")),
              expression("m"^"6"*"A Site"),
              expression("Literature")),
      cex = 0.8,
      col = 'black'
    ),
    points = list(
      pch = c(20, 15, 18),
      col = c('transparent', 'black', 'black'),
      cex = c(0, 1, 1),
      alpha = 1
    ),
    x = 0.85,
    y = 0.95,
    padding.text = 1
  ),
  description = 'Manhattan plot created using BoutrosLab.plotting.general'
)
system(paste0("cp ", filename, " ~/figures"))
