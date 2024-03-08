source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")

library(RNOmni)

# Preamble ----------------------------------------------------------------
# Makes boxplots for Protein
# Saves Protein linear models

# Loading Data ------------------------------------------------------------

# Setting variables
meth = "MeTPeak"
tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"
output.tag = 'CPCGENE_filtered.final'
n_components = '0'

# Sinha et al., 2019 - Protein Data
protein = read.csv("Summary_Tables/Sinha.2019.data/Sinha_2019_Protein.csv",
                   header = T, skip = 2, stringsAsFactors = F)

# Loading Genotype Data
# geno.file = paste0("Input/", output.tag, ".matrix.tsv")
# geno = read.delim(geno.file, header = T)

# Loading Covariate Data
cov.file = "ComplementaryQTL/protein.covariates"
cov.table = read.delim(cov.file, header = T)

# Significant results
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.rsav")))
# [1] "sig.results"

# Protein samples
protein.samples = readLines("ComplementaryQTL/protein_sample_file.tsv")

# Analysis ----------------------------------------------------------------

changeSciNot <- function(n) {
  output <- formatC(n, digits = 2, format = "e")
  num = strsplit(output, split = "e")[[1]][1]
  exponent = strtoi(gsub("-0", "-", strsplit(output, split = "e")[[1]][2]))
  if(num == "1.00"){
    output = bquote("P =" ~ 10^.(exponent))
  }else{
    output = bquote("P =" ~ .(num) ~ x ~ 10^.(exponent))
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

  # Going from peak to gene
  gene = gene_conversion[gsub("[:].*", "", gene)]

  # If there's no detected protein
  if(!gene %in% protein$Gene){
    return(data.frame("protein.pval" = NA, "protein.effectSize" = NA))
  }

  # Extracting info
  geno.snp = as.numeric(geno[geno$id == snp, protein.samples])
  expr.snp = protein[protein$Gene == gene, protein.samples]

  # If there are multiple proteins groups for that gene,
  # Take the one with the fewest NA's or the greatest mean abundance
  if(nrow(expr.snp) > 1){
    check.na = apply(expr.snp, 1, function(x) sum(is.na(x)))
    expr.snp = expr.snp[which.min(check.na),]
    if(nrow(expr.snp) > 1){
      check.mean = apply(expr.snp, 1, mean)
      expr.snp = expr.snp[which.max(check.mean),]
    }
  }
  expr.snp = as.numeric(expr.snp)

  # Non NA's
  keep = which(!is.na(expr.snp))
  geno.snp = geno.snp[keep]
  expr.snp = RNOmni::RankNorm(expr.snp[keep])
  samples = protein.samples[keep]

  # Creating Plotting Data
  plt.df = data.frame(
    # "Sample" = samples,
    "Genotype" = geno.snp,
    "Expression" = expr.snp,
    t(cov.table[,samples]),
    stringsAsFactors = F
  )
  # Recreating linear model
  if(corrected){
    corrected.lm = lm(formula = Expression ~ ., data = plt.df)
    corrected.expr = plt.df$Expression - rowSums(coef(corrected.lm)[-(1:2)]*plt.df[,3:ncol(plt.df)])
    plt.df$Expression = corrected.expr
    beta.val = coef(corrected.lm)["Genotype"]
    p.val = summary(corrected.lm)$coefficients["Genotype", "Pr(>|t|)"]
  } else {
    uncorrected.lm = lm(formula = Expression ~ Genotype, data = plt.df)
    beta.val = coef(uncorrected.lm)["Genotype"]
    p.val = summary(uncorrected.lm)$coefficients["Genotype", "Pr(>|t|)"]
  }

  # Converting info
  geno.conv = structure(c("AA", "AB", "BB"), names = as.character(c(0, 1, 2)))
  plt.df$Genotype = geno.conv[as.character(plt.df$Genotype)]
  plt.df$Genotype = factor(plt.df$Genotype, levels = c("AA", "AB", "BB"))

  # Creating Title
  snp.main = gsub("[:].*", "", snp)
  combined.main = paste0(snp.main, " | ", gene)

  # Legend
  legend.stats <- legend.grob(
    list(
      legend = list(
        colours = "transparent",
        title = changeSciNot(p.val),
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
    ylab.label = expression(bold(paste("Protein"))),
    xlab.label = "Genotype",
    ylab.cex = 1.5,
    xlab.cex = 1.5,
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xaxis.lab = xaxis_lab,
    # ylab.axis.padding = 1.5,
    yaxis.tck = 0,
    xaxis.tck = 0,
    right.padding = 5,
    # Legend
    legend = list(
      inside = list(
        x = 0.5,
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

  return(data.frame("protein.pval" = p.val, "protein.effectSize" = beta.val))
}

# Plotting ----------------------------------------------------------------

filename = paste0("~/figures/", meth, ".", tag, ".", output.tag, "_",  n_components, ".protein.corrected.pdf")
pdf(filename, width = 5, height = 6)
protein.results = data.frame()
for(i in 1:nrow(sig.results)){
  tmp = create.plots(
    snp = sig.results$snps[i],
    gene =  sig.results$gene[i],
    corrected = T)
  protein.results = rbind(protein.results, tmp)
}
dev.off()

protein.results = cbind.data.frame(sig.results[,c("snps", "gene", "GeneName", "Literature")], protein.results)
protein.results$protein.fdr = p.adjust(protein.results$protein.pval, method = "fdr")
protein.results = protein.results[order(protein.results$protein.pval),]

# Saving Protein Results --------------------------------------------------

save(protein.results,
     file = paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.protein.rsav"))
