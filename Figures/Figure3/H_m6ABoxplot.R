source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")

library("MatrixEQTL")

# Preamble ----------------------------------------------------------------
# Generating Individual Boxplots Showing the Association

# Loading Data ------------------------------------------------------------

# Setting variables
meth = "MeTPeak"
tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"
output.tag = 'CPCGENE_filtered.final'
n_components = '0'

# Loading Expression Data
expr.file = paste0("Input/", meth, ".", tag, ".Expression.tsv")
expr = read.delim(expr.file, header = T)

# Loading Genotype Data
# geno.file = paste0("Input/", output.tag, ".matrix.tsv")
# geno = read.delim(geno.file, header = T)

# Loading Covariate Data
cov.file = paste0("covariates/", meth, ".", tag, "_", n_components, "_peer.tsv")
cov.table = read.delim(cov.file, header = T)

# Significant results
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.rsav")))
# [1] "sig.results" "geno"

# Making Plotting Matrix --------------------------------------------------

# changeSciNot <- function(p, b) {
#   output <- formatC(p, digits = 2, format = "e")
#   num = strsplit(output, split = "e")[[1]][1]
#   exponent = strtoi(gsub("-0", "-", strsplit(output, split = "e")[[1]][2]))
#   effectsize = formatC(b, digits = 3)
#   if(num == "1.00"){
#     output = bquote(atop("Q =" ~ 10^.(exponent), .(sym("beta")) ~ " = " ~ .(effectsize)))
#   }else{
#     output = bquote(atop("Q =" ~ .(num) ~ x ~ 10^.(exponent), .(sym("beta")) ~ " = " ~ .(effectsize)))
#   }
#   output
# }

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

create.plots = function(snp, gene, corrected = T){

  # snp = sig.results$snps[1]
  # gene = sig.results$gene[1]
  cat(snp, ", ", gene, "\n")

  # xaxis labels
  snp_split <- strsplit(snp, split = ":")[[1]]
  xaxis_lab <- c(
    paste0(snp_split[4], snp_split[4]),
    paste0(snp_split[4], snp_split[5]),
    paste0(snp_split[5], snp_split[5])
  )

  # Extracting info
  geno.snp = as.numeric(geno[geno$id == snp, all.euro.samples])
  expr.snp = as.numeric(expr[expr$peak == gene, all.euro.samples])
  beta.val = sig.results$beta[sig.results$gene == gene & sig.results$snps == snp]
  adjp.val = sig.results$FDR[sig.results$gene == gene & sig.results$snps == snp]

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
  gene.main = gene_conversion[gsub(":.*", "", gene)]
  peak.main = gsub(".*:", "", gene)
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

filename = paste0("~/figures/", meth, ".", tag, ".", output.tag, "_",  n_components, ".corrected.pdf")
pdf(filename, width = 5, height = 6)
for(i in 1:nrow(sig.results)){
  create.plots(
    snp = sig.results$snps[i],
    gene =  sig.results$gene[i],
    corrected = T)
}
dev.off()
