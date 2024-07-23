source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")

library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.survival);

# Preamble ----------------------------------------------------------------
# This tests for clinical associations with SNPs identified to be QTLs

# Loading Data ------------------------------------------------------------

# Setting variables
meth = "MeTPeak"
tag = "V9_MetricVoting_OptimizationFunction_MaxGapOnly"
output.tag = 'CPCGENE_filtered.final'
n_components = '0'

# Loading Covariate Data
cov.file = paste0("covariates/", meth, ".", tag, "_", n_components, "_peer.tsv")
cov.table = read.delim(cov.file, header = T)

# Significant results
print(load(paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.rsav")))
# [1] "sig.results" "geno"

# Clinical Data -----------------------------------------------------------

# Clinical Data
clin.vars = c("pathologic_isup_grade", "summary_t", "pre_treatment_psa")
filename = "Summary_Tables/2019-08-20_500pg_SupTable1_feature_by_patient.tsv"
clinical = read.delim(filename, header = T)
rownames(clinical) = clinical$sample
clinical = clinical[all.euro.samples, clin.vars]

# Grouping ISUP
clinical$pathologic_isup_grade = as.character(clinical$pathologic_isup_grade)
clinical$pathologic_isup_grade[clinical$pathologic_isup_grade %in% c("4", "5")] <- "4+"

table(clinical$pathologic_isup_grade)
# 1  2  3 4+
#   9 82 30 11

# Grouping Clinical T Category
clinical$summary_t = substr(clinical$summary_t, 1, 2)

table(clinical$summary_t)
# T2 T3
# 70 62

# IDC status
filename = "Summary_Tables/2021-05-07_m6A_idc_or_cribiform.txt"
idc = read.delim(filename, header = T)
clinical$idc = as.numeric(idc$idc_or_cribiform[match(rownames(clinical), idc$patient_id)])

table(clinical$idc)
# 0  1
# 70 48

# Making Genotypes Binary -------------------------------------------------

# Reformatting
geno.df = geno
rownames(geno.df) = geno.df$id
geno.df = geno.df[,all.euro.samples]

# Making it binary
geno.df[geno.df == 1] <- 0
geno.df[geno.df == 2] <- 1

# Analysis ----------------------------------------------------------------

format.results = function(df, tag){
  df = df[,match(colnames(df), c('p', 'loci', 'var'))]
  df$fdr = p.adjust(df$p, method = 'fdr')
  df = df[,c("loci", "p", 'fdr')]
  colnames(df)[colnames(df) %in% c('p', 'fdr')] <- paste0(tag, ".", colnames(df)[colnames(df) %in% c('p', 'fdr')])
  df[order(df[,paste0(tag, ".p")]),]
}

# Categorical Variables
categ = function(df, vari){

  # Check column names
  all(colnames(df) == rownames(clinical))

  # Clinical variable
  clin.vec = clinical[,vari]
  sample.keep = !is.na(clin.vec)
  clin.vec = clin.vec[sample.keep]

  # Filtering rows that only have one value
  df = df[,sample.keep]
  keep = apply(df, 1, function(x) length(unique(x[!is.na(x)])) > 1)
  df = df[keep,]

  # Initializing results
  df = t(df)
  res = as.data.frame(
    matrix(
      nrow = ncol(df),
      ncol = 3,
      byrow = TRUE
    ),
    stringsAsFactors = FALSE
  );

  colnames(res) <- c('p', 'loci', 'var');
  rownames(res) <- colnames(df);

  # Chi-squared results
  for(i in colnames(df)){
    cat(i, "\n")
    test.res = chisq.test(
      x = as.character(df[,i]),
      y = clin.vec
    )
    res[i, "loci"] <- i
    res[i, "var"] <- vari
    res[i, "p"] <- test.res$p.value
  }
  return(res)
}

# ISUP Results
ISUP.results = categ(df = geno.df, vari = "pathologic_isup_grade")
ISUP.results = format.results(ISUP.results, "ISUP")

# Summary T Results
Summary_T.results = categ(df = geno.df, vari = "summary_t")
Summary_T.results = format.results(Summary_T.results, "T")

# IDC Results
IDC.results = categ(df = geno.df, vari = "idc")
IDC.results = format.results(IDC.results, "IDC")

# Continuous Variables
cont = function(df, vari){

  # Check column names
  all(colnames(df) == rownames(clinical))

  # Clinical variable
  clin.vec = clinical[,vari]
  sample.keep = !is.na(clin.vec)
  clin.vec = clin.vec[sample.keep]

  # Filtering rows that only have one value
  df = df[,sample.keep]
  keep = apply(df, 1, function(x) length(unique(x[!is.na(x)])) > 1)
  df = df[keep,]

  df = t(df)
  res = as.data.frame(
    matrix(
      nrow = ncol(df),
      ncol = 3,
      byrow = TRUE
    ),
    stringsAsFactors = FALSE
  );

  colnames(res) <- c('p', 'loci', 'var');
  rownames(res) <- colnames(df);

  for(i in colnames(df)){
    cat(i, "\n")
    x = clinical[,vari][df[,i] == 1]
    y = clinical[,vari][df[,i] == 0]
    test.res = wilcox.test(
      x = x,
      y = y
    )
    res[i, "loci"] <- i
    res[i, "var"] <- vari
    res[i, "p"] <- test.res$p.value
  }
  return(res)
}

# PSA
PSA.results = cont(geno.df, "pre_treatment_psa")
PSA.results = format.results(PSA.results, "PSA")

clinical.results = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "loci", all = TRUE),
       list(ISUP.results, Summary_T.results, IDC.results, PSA.results))

# Saving Clinical Results -------------------------------------------------

save(clinical.results,
     file = paste0("results/", meth, ".",  tag, ".", output.tag, ".", n_components, ".MatrixQTL.cis.sig.clinical.rsav"))

# Making Some Contingency Tables ------------------------------------------

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

make.contigency.table = function(df, this.main, rsid){
  at.cols = seq(0, max(c(rowSums(df), colSums(df)))+10, 10)
  middle.heatmap = create.heatmap(
    df,
    clustering.method = "none",
    same.as.matrix = T,
    # Colourkey
    colour.scheme = c("white", "dodgerblue"),
    at = at.cols,
    print.colour.key = F,
    # Axes
    axes.lwd = 1,
    ylab.label = expression("ISUP Grade Group"),
    ylab.cex = 1,
    yaxis.cex = 1,
    yaxis.lab = rownames(df),
    yaxis.fontface = 1,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.rot = 0,
    # Adding Text
    cell.text = as.vector(df),
    row.pos = rep(nrow(df):1, ncol(df)),
    col.pos = sort(rep(1:ncol(df), nrow(df))),
    text.fontface = 1,
    text.cex = 1,
    text.col = 'black',
    text.position = NULL,
    text.offset = 0
  )

  right.hm.data = data.frame("x" = rowSums(df), "y" = rowSums(df))
  right.heatmap = create.heatmap(
    right.hm.data,
    clustering.method = "none",
    same.as.matrix = T,
    # Colourkey
    colour.scheme = c("white", "dodgerblue"),
    at = at.cols,
    print.colour.key = F,
    # Axes
    axes.lwd = 1,
    xaxis.tck = 0,
    yaxis.tck = 0,
    # Adding Text
    row.pos = nrow(right.hm.data):1,
    col.pos = rep(1.5, nrow(right.hm.data)),
    cell.text = right.hm.data$x,
    text.fontface = 1,
    text.cex = 1,
    text.col = 'black',
    text.position = NULL,
    text.offset = 0
  )

  bottom.hm.data = data.frame("x" = colSums(df), "y" = colSums(df))
  bottom.heatmap = create.heatmap(
    bottom.hm.data,
    clustering.method = "none",
    same.as.matrix = F,
    # Colourkey
    colour.scheme = c("white", "dodgerblue"),
    at = at.cols,
    print.colour.key = F,
    # Axes
    axes.lwd = 1,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.lab = c("AA + AB", "BB"),
    xaxis.cex = 1,
    xlab.label =  as.expression(bquote(.(rsid))),
    xlab.cex = 1,
    xaxis.rot = 0,
    xaxis.fontface = 1,
    # Adding Text
    row.pos = rep(1.5, nrow(bottom.hm.data)),
    col.pos = 1:nrow(bottom.hm.data),
    cell.text = bottom.hm.data$x,
    text.fontface = 1,
    text.cex = 1,
    text.col = 'black',
    text.position = NULL,
    text.offset = 0
  )

  plt = create.multipanelplot(
    plot.objects = list(middle.heatmap, right.heatmap, bottom.heatmap),
    layout.width = 2,
    layout.height = 2,
    main = this.main,
    main.cex = 1,
    main.x = 0.5,
    main.y = 0,
    x.spacing = -2,
    y.spacing = -5.5,
    ylab.axis.padding = -1,
    xlab.axis.padding = 1,
    layout.skip = c(F, F, F, T),
    plot.objects.heights = c(5, 2),
    plot.objects.widths = c(5, 1.5)
  )
  return(plt)
}

loci = c("rs143089027:chr1:162583262:C:T", "rs76338659:chr13:61422806:T:C", "rs57557217:chr5:116581824:A:G")


filename = "~/figures/340_Clinical.ISUP.Geno.pdf"
pdf(filename, width = 4, height = 5)
for(i in loci){
  df = data.frame(
    "ISUP" = clinical[all.euro.samples, "pathologic_isup_grade"],
    "Geno" = as.numeric(geno.df[i, all.euro.samples])
  )
  df = table(df$ISUP, df$Geno)
  fdr = ISUP.results$ISUP.fdr[ISUP.results$loci == i]
  print(make.contigency.table(df, this.main = changeSciNot(fdr), rsid = gsub("[:].*", "", i)))
}

dev.off()


# Checking whether SNPs are ancestry related ------------------------------

make.ancestry.contigency.table = function(df, rsid_1, rsid_2){
  at.cols = seq(0, max(c(rowSums(df), colSums(df)))+10, 10)
  middle.heatmap = create.heatmap(
    df,
    clustering.method = "none",
    same.as.matrix = T,
    # Colourkey
    colour.scheme = c("white", "dodgerblue"),
    at = at.cols,
    print.colour.key = F,
    # Axes
    axes.lwd = 1,
    ylab.label = as.expression(bquote(.(rsid_1))),
    ylab.cex = 1,
    yaxis.cex = 1,
    yaxis.lab = c("AA", "AB", "BB"),
    yaxis.fontface = 1,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.rot = 0,
    # Adding Text
    cell.text = as.vector(df),
    row.pos = rep(nrow(df):1, ncol(df)),
    col.pos = sort(rep(1:ncol(df), nrow(df))),
    text.fontface = 1,
    text.cex = 1,
    text.col = 'black',
    text.position = NULL,
    text.offset = 0
  )

  right.hm.data = data.frame("x" = rowSums(df), "y" = rowSums(df))
  right.heatmap = create.heatmap(
    right.hm.data,
    clustering.method = "none",
    same.as.matrix = T,
    # Colourkey
    colour.scheme = c("white", "dodgerblue"),
    at = at.cols,
    print.colour.key = F,
    # Axes
    axes.lwd = 1,
    xaxis.tck = 0,
    yaxis.tck = 0,
    # Adding Text
    row.pos = nrow(right.hm.data):1,
    col.pos = rep(1.5, nrow(right.hm.data)),
    cell.text = right.hm.data$x,
    text.fontface = 1,
    text.cex = 1,
    text.col = 'black',
    text.position = NULL,
    text.offset = 0
  )

  bottom.hm.data = data.frame("x" = colSums(df), "y" = colSums(df))
  bottom.heatmap = create.heatmap(
    bottom.hm.data,
    clustering.method = "none",
    same.as.matrix = F,
    # Colourkey
    colour.scheme = c("white", "dodgerblue"),
    at = at.cols,
    print.colour.key = F,
    # Axes
    axes.lwd = 1,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.lab = c("AA", "AB", "BB"),
    xaxis.cex = 1,
    xlab.label = as.expression(bquote(.(rsid_2))),
    xlab.cex = 1,
    xaxis.rot = 0,
    xaxis.fontface = 1,
    # Adding Text
    row.pos = rep(1.5, nrow(bottom.hm.data)),
    col.pos = 1:nrow(bottom.hm.data),
    cell.text = bottom.hm.data$x,
    text.fontface = 1,
    text.cex = 1,
    text.col = 'black',
    text.position = NULL,
    text.offset = 0
  )

  plt = create.multipanelplot(
    plot.objects = list(middle.heatmap, right.heatmap, bottom.heatmap),
    layout.width = 2,
    layout.height = 2,
    # main = this.main,
    # main.cex = 1,
    # main.x = 0.5,
    # main.y = 0,
    x.spacing = -3,
    y.spacing = -5.5,
    ylab.axis.padding = 1,
    xlab.axis.padding = 1,
    layout.skip = c(F, F, F, T),
    plot.objects.heights = c(5, 2.5),
    plot.objects.widths = c(5, 1.5)
  )
  return(plt)
}

# Making some contingency tables
geno_ancestry = geno
rownames(geno_ancestry) = geno_ancestry$id
geno_ancestry = geno_ancestry[loci, all.euro.samples]
rownames(geno_ancestry) = gsub("[:].*", "", rownames(geno_ancestry))

plt1 = make.ancestry.contigency.table(
  df = table(
    geno_ancestry["rs143089027",],
    geno_ancestry["rs76338659",]
  ),
  rsid_1 = "rs143089027",
  rsid_2 = "rs76338659"
)
# 0  1  2
# 0 77 35  2
# 1  8  8  0
# 2  1  1  0

plt2 = make.ancestry.contigency.table(
  df = table(
    geno_ancestry["rs143089027",],
    geno_ancestry["rs57557217",]
  ),
  rsid_1 = "rs143089027",
  rsid_2 = "rs57557217"
)
# 0  1  2
# 0 60 47  7
# 1  8  8  0
# 2  1  0  1

plt3 = make.ancestry.contigency.table(
  df = table(
    geno_ancestry["rs76338659",],
    geno_ancestry["rs57557217",]
  ),
  rsid_1 = "rs76338659",
  rsid_2 = "rs57557217"
)
# 0  1  2
# 0 41 39  6
# 1 26 16  2
# 2  2  0  0

filename = "~/figures/340_Clinical.Ancestry.pdf"
pdf(filename, width = 4, height = 4)

print(plt1)
print(plt2)
print(plt3)

dev.off()
