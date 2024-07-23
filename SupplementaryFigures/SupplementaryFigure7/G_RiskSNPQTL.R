source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")



# Preamble ----------------------------------------------------------------
# Checking local QTLs for risk SNPs
# Inspired by Yuan et al., 2022

# Loading Data ------------------------------------------------------------

# Risk SNPs
risk_snps = read.delim("Input/CPCGENE_RiskSNPs.final.location.tsv", header = T)

# Local QTL data
print(load("results/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.CPCGENE_filtered.final.0.all.rsav"))
# [1] "me"

# Cis eqtls
cis.qtls = me$cis$eqtls

# Loading Expression Data
expr.file = "Input/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.Expression.tsv"
expr = read.delim(expr.file, header = T)

# Loading Genotype Data
geno.file = "Input/CPCGENE_RiskSNPs.final.matrix.tsv"
geno = read.delim(geno.file, header = T)

# Loading Covariate Data
cov.file = "covariates/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly_0_peer.tsv"
cov.table = read.delim(cov.file, header = T)

# Analysis ----------------------------------------------------------------

# Extracting risk SNP QTL
risk_snp_cis = cis.qtls[cis.qtls$snps %in% risk_snps$ID,]

# Recomputing FDR
risk_snp_cis$FDR = p.adjust(risk_snp_cis$pvalue, method = "fdr")

dim(risk_snp_cis[risk_snp_cis$FDR < 0.1,])
# [1] 13  6

# Genes
colnames(risk_snp_cis)[colnames(risk_snp_cis) == "gene"] <- "peak"
risk_snp_cis$gene = gsub("[:].*", "", risk_snp_cis$peak)
risk_snp_cis$gene_name = gene_conversion[risk_snp_cis$gene]

# Generating a set of significant results
sig.results = risk_snp_cis[risk_snp_cis$FDR < 0.1,]

# Generating significant results stable -----------------------------------

write.table(
  sig.results,
  file = "~/figures/local_qtl_significant_stable.tsv",
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)


# Summary Plots -----------------------------------------------------------

# This is currently manually chosen based on what's in the paper
# Should modify if ever need
label.results <- sig.results[c(1, 4),]
label.results$rsid <- gsub("[:].*", "", label.results$snps)
label.results$peak_id <- gsub(".*[:]", "", label.results$peak)
label.results$label <- paste0(label.results$rsid, " | ", label.results$gene_name, " | ", label.results$peak_id)

filename = "~/figures/511_local_risk_snps_QTL.pdf"
pdf(filename, width = 5, height = 5)

# Histogram of P-values
create.histogram(
  risk_snp_cis$pvalue,
  breaks = seq(0, 1, 0.1),
  # type = "count",
  # Axes
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xlab.cex = 1.5,
  ylab.cex = 1.5,
  xlab.label = expression(bold("Local risk SNP m"^"6"*"A QTL P value")),
  ylab.label = "Percent",
  xaxis.tck = 0,
  yaxis.tck = 0
  # main.cex = 1.5,
  # main = expression(bold("Local risk SNP m"^"6"*"A QTL"))
)

# Volcano plot
create.scatterplot(
  -log10(pvalue) ~ beta,
  risk_snp_cis,
  # Axes
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xlab.cex = 1.5,
  ylab.cex = 1.5,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xlab.label = expression(bold("Local risk SNP m"^"6"*"A QTL "*beta)),
  ylab.label = expression(bold("-log"[10]*" P value")),
  xlimits = c(-1, 1),
  ylimits = c(0, 5.5),
  xat = seq(-1, 1, 1),
  yat = seq(0, 5, 1),
  # main = expression(bold("Local risk SNP m"^"6"*"A QTL")),
  # main.cex = 1.5,
  # Adding sig points
  add.points = TRUE,
  points.x = label.results$beta,
  points.y = -log10(label.results$pvalue),
  points.pch = 19,
  points.col = 'blue',
  points.col.border = 'blue',
  points.cex = 1,
  # Adding labels
  add.text = TRUE,
  text.labels = label.results$label,
  text.x = label.results$beta + 0.5,
  text.y = -log10(label.results$pvalue) + 0.2,
  text.col = 'blue',
  text.cex = 1.2,
  text.fontface = 'bold',
  # text.guess.labels = FALSE,
  # text.guess.skip.labels = TRUE,
  # text.guess.ignore.radius = FALSE,
  # text.guess.ignore.rectangle = FALSE,
  # text.guess.radius.factor = 1,
  # text.guess.buffer.factor = 1,
  # text.guess.label.position = NULL,
)

dev.off()

# Individual Plots --------------------------------------------------------

# Boxplot of the QTLs
changeSciNot <- function(n) {
  output <- formatC(n, digits = 2, format = "e")
  num = strsplit(output, split = "e")[[1]][1]
  exponent = strtoi(gsub("-0", "-", strsplit(output, split = "e")[[1]][2]))
  if(num == "1.00"){
    output = bquote("Q =" ~ 10^.(exponent))
  }else{
    output = bquote("Q =" ~ .(num) ~ x ~ 10^.(exponent))
  }
  output
}

create.plots = function(snp, peak, corrected = T){

  # snp = sig.results$snps[1]
  # peak = sig.results$peak[1]
  cat(snp, ", ", peak, "\n")

  # xaxis labels
  snp_split <- strsplit(snp, split = ":")[[1]]
  xaxis_lab <- c(
    paste0(snp_split[4], snp_split[4]),
    paste0(snp_split[4], snp_split[5]),
    paste0(snp_split[5], snp_split[5])
  )

  # Extracting info
  geno.snp = as.numeric(geno[geno$id == snp, all.euro.samples])
  expr.snp = as.numeric(expr[expr$peak == peak, all.euro.samples])
  beta.val = sig.results$beta[sig.results$peak == peak & sig.results$snps == snp]
  adjp.val = sig.results$FDR[sig.results$peak == peak & sig.results$snps == snp]

  # Creating Plotting Data
  plt.df = data.frame(
    # "Sample" = colnames(geno.snp),
    "Genotype" = geno.snp,
    "Expression" = expr.snp,
    t(cov.table),
    stringsAsFactors = F
  )
  # Recreating linear model
  if(corrected){
    corrected.lm = lm(formula = Expression ~ ., data = plt.df)
    corrected.expr = plt.df$Expression - rowSums(coef(corrected.lm)[-(1:2)]*plt.df[,3:ncol(plt.df)])
    plt.df$Expression = corrected.expr
  } else {
    uncorrected.lm = lm(formula = Expression ~ Genotype, data = plt.df)
  }

  # Converting info
  geno.conv = structure(c("AA", "AB", "BB"), names = as.character(c(0, 1, 2)))
  plt.df$Genotype = geno.conv[as.character(plt.df$Genotype)]
  plt.df$Genotype = factor(plt.df$Genotype, levels = c("AA", "AB", "BB"))

  # Creating Title
  snp.main = gsub("[:].*", "", snp)
  gene.main = gene_conversion[gsub(":.*", "", peak)]
  peak.main = gsub(".*:", "", peak)
  combined.main = paste0(snp.main, " | ", gene.main, " | ", peak.main)

  # Legend
  legend.stats <- legend.grob(
    list(
      legend = list(
        colours = "transparent",
        title = changeSciNot(adjp.val),
        labels = "",
        border = 'transparent'
      )
    ),
    size = 1,
    title.cex = 1.5,
    title.just = 'left'
  );


  bp = create.boxplot(
    Expression ~ Genotype,
    plt.df,
    # filename = "~/figures/test2.png",
    # Main
    main = combined.main,
    main.cex = 1.5,
    main.just = "center",
    # main.x = 0.25,
    # main.y = 0.5,
    # Strip plot
    add.stripplot = TRUE,
    jitter.factor = 1,
    points.pch = 19,
    points.col = "darkgrey",
    points.cex = 1,
    points.alpha = 0.5,
    # Axes
    ylab.label = expression(bold(paste("m"^"6"*"A"))),
    xlab.label = "Genotype",
    ylab.cex = 1.5,
    xlab.cex = 1.5,
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xaxis.lab = xaxis_lab,
    ylimits = c(min(plt.df$Expression)-1, max(plt.df$Expression) + 1.5),
    # ylab.axis.padding =0.8,
    yaxis.tck = 0,
    xaxis.tck = 0,
    right.padding = 5,
    # Legend
    legend = list(
      inside = list(
        x = 0.05,
        y = 0.97,
        size = 50,
        fun = legend.stats
      )
    )# ,
    # plot style
    # height = 6,
    # width = 5,
    # resolution = 1600
  )

  print(bp)

}

# Plotting ----------------------------------------------------------------

filename = "~/figures/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.CPCGENE_RiskSNPs.0.corrected.pdf"

pdf(filename, width = 5, height = 6)
for(i in 1:nrow(sig.results)){
  create.plots(
    snp = sig.results$snps[i],
    peak =  sig.results$peak[i],
    corrected = T)
}
dev.off()

filename = "results/riskSNP_local_QTL.tsv"
write.table(
  sig.results,
  filename,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
system(paste0("cp ", filename, " ~/figures"))
