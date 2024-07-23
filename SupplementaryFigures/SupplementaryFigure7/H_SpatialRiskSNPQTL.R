source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")


# Preamble ----------------------------------------------------------------
# Makes boxplots for RNA for sig risk SNP QTLs
# Saves RNA linear models

# Loading Data ------------------------------------------------------------

# Chen et al., 2019 - RNAseq data
bulk.file = "ComplementaryQTL/bulkRNA.Expression.tsv"
bulk.RNA = read.delim(bulk.file, header = T)

# Loading Genotype Data
geno.file <- "Input/CPCGENE_RiskSNPs.final.matrix.tsv"
geno <- read.delim(geno.file, header = T)

# Loading Covariate Data
cov.file = "ComplementaryQTL/bulkRNA.covariates"
cov.table = read.delim(cov.file, header = T)

# Significant results
sig.results <- read.delim("results/riskSNP_local_QTL.tsv", header = T)

# RNA samples
rna.samples = readLines("ComplementaryQTL/bulkRNA_sample_file.tsv")

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
  # snp = "rs4876:chr11:32103432:T:C"
  # gene = "ENSG00000285283.1"
  cat(snp, ", ", gene, "\n")

  # xaxis labels
  snp_split <- strsplit(snp, split = ":")[[1]]
  xaxis_lab <- c(
    paste0(snp_split[4], snp_split[4]),
    paste0(snp_split[4], snp_split[5]),
    paste0(snp_split[5], snp_split[5])
  )

  # Going from peak to gene
  gene = gsub("[:].*", "", gene)

  # Extracting info
  geno.snp = as.numeric(geno[geno$id == snp, rna.samples])
  expr.snp = as.numeric(bulk.RNA[bulk.RNA$gene_id == gene, rna.samples])

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
  gene.main = gene_conversion[gsub(":.*", "", gene)]
  combined.main = paste0(snp.main, " | ", gene.main)

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
    ylab.label = expression(bold(paste("RNA"))),
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

  return(data.frame("bulkRNA.pval" = p.val, "bulkRNA.effectSize" = beta.val))
}

# Plotting ----------------------------------------------------------------

filename = "~/figures/riskSNP_local_QTL.RNA.corrected.pdf"
pdf(filename, width = 5, height = 6)
rna.results = data.frame()
for(i in 1:nrow(sig.results)){
  tmp = create.plots(
    snp = sig.results$snps[i],
    gene =  sig.results$gene[i],
    corrected = T)
  rna.results = rbind(rna.results, tmp)
}
dev.off()

rna.results = cbind.data.frame(sig.results[,c("snps", "gene", "gene_name")], rna.results)
rna.results$bulkRNA.fdr = p.adjust(rna.results$bulkRNA.pval, method = "fdr")
rna.results = rna.results[order(rna.results$bulkRNA.pval),]

# Saving RNA Results ------------------------------------------------------

save(rna.results, file = "results/riskSNP_local_QTL.MatrixQTL.cis.sig.RNA.rsav")

filename <- "~/figures/riskSNP_local_QTL.RNA.tsv"
write.table(
  rna.results,
  file = filename,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
