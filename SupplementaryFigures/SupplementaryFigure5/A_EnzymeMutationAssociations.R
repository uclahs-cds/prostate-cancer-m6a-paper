source("/cluster/home/helenzhu/code/snakemake_M6A_34_EnzymeMutations/000_HEADER.R")

library(BoutrosLab.plotting.general)


# Preamble ----------------------------------------------------------------
# This makes a dotmap of the mutation associations
# RNA
# Protein
# Number of Peaks
# Clinical Variables

# Loading Data ------------------------------------------------------------

# RNA
print(load("M_M6A_Enzyme_Mutations/m6a.enzymes.rna.diff.data.rsav"))
# [1] "diff.rna"  "chen"      "mutations"

# Protein
print(load("M_M6A_Enzyme_Mutations/m6a.enzymes.protein.diff.data.rsav"))
# [1] "differential.prot" "protein"           "mutations"

# Number of Peaks
print(load("M_M6A_Enzyme_Mutations/m6a.enzymes.peaks.diff.data.rsav"))
# [1] "differential.peaks" "metpeak"            "mutations"

# Clinical Variables
print(load("M_M6A_Enzyme_Mutations/m6a.enzymes.clinical.data.rsav"))
# [1] "surv.results"      "ISUP.results"      "IDC.results"
# [4] "Summary_T.results" "Age.results"       "PSA.results"
# [7] "Hypoxia.results"   "bcr_data"          "clinical"
# [10] "mutations"

# Making an Overview Plot -------------------------------------------------

# RNA
diff.rna = diff.rna[,c("gene", "cna", "class", "pval", "fc", "fdr")]
colnames(diff.rna) = c("Gene", "Mutation", "Class", "RNA.p", "RNA.effectSize", "RNA.fdr")

# Protein
differential.prot = differential.prot[,c("gene", "cna", "class", "pval", "fc", "fdr")]
colnames(differential.prot) = c("Gene", "Mutation", "Class", "Protein.p", "Protein.effectSize", "Protein.fdr")

# Peaks
differential.peaks = differential.peaks[,c("gene", "cna", "class", "metpeak.pval", "metpeak.fc", "metpeak.fdr")]
colnames(differential.peaks) = c("Gene", "Mutation", "Class", "MeTPeak.p", "MeTPeak.effectSize", "MeTPeak.fdr")

# Do something about surv.results
colnames(surv.results) = c("Gene", "Mutation", "Class", "BCR.effectSize", "BCR.p", "BCR.fdr")

# Compiling Plotting Data
df.list = list(
  diff.rna,
  differential.prot,
  differential.peaks,
  ISUP.results,
  IDC.results,
  Summary_T.results,
  Age.results,
  PSA.results,
  Hypoxia.results,
  surv.results)

df.compiled = Reduce(function(df1, df2) merge(df1, df2, by = c("Gene", "Mutation", "Class"), all = TRUE), df.list)

# Write an output file for Rupert
write.table(
  df.compiled,
  file = "~/m6A_enzyme_mutations.tsv",
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

# P value
df.Pvalue = df.compiled[,grep(".p", colnames(df.compiled), value = T, fixed = T)]
colnames(df.Pvalue) = gsub(".p", "", colnames(df.Pvalue), fixed = T)

# Q value
df.Qvalue = df.compiled[,grep(".fdr", colnames(df.compiled), value = T, fixed = T)]
colnames(df.Qvalue) = gsub(".fdr", "", colnames(df.Qvalue), fixed = T)

# Fold Change
df.foldChange = df.compiled[,grep(".effectSize", colnames(df.compiled), value = T, fixed = T)]
df.foldChange = df.foldChange[,!colnames(df.foldChange) %in% c("Age.effectSize", "PSA.effectSize")] # Currently, this isn't a meaningful distinction
df.foldChange[,c("RNA.effectSize", "Protein.effectSize", "MeTPeak.effectSize", "BCR.effectSize")] = log2(df.foldChange[,c("RNA.effectSize", "Protein.effectSize", "MeTPeak.effectSize", "BCR.effectSize")])
colnames(df.foldChange) = gsub(".effectSize", "", colnames(df.foldChange), fixed = T)

missing.cols = setdiff(colnames(df.Qvalue), colnames(df.foldChange))
df.fc.missing = matrix(0, ncol = length(missing.cols), nrow = nrow(df.foldChange), dimnames = list(rownames(df.foldChange), missing.cols))
df.foldChange = cbind(df.foldChange, df.fc.missing)
df.foldChange = df.foldChange[,colnames(df.Qvalue)]

# Check
all(colnames(df.foldChange) == colnames(df.Pvalue))
all(rownames(df.foldChange) == rownames(df.Pvalue))

# Row names
rnames = df.compiled$gene_name

# Col names
cnames = colnames(df.Pvalue)

# Creating a Dotmap -------------------------------------------------------

# Fold Change
col.foldChange.func = function(i){
  cols = rep("darkorange1", length(i))
  cols[i < 0] <- "dodgerblue2"
  cols[i == 0] <- "transparent"
  return(cols)
}

key.sizes = c(-1.0, -0.5, -0.1, 0, 0.1, 0.5, 1.0)
key.foldChange.func = function(x) { 0.1 + (2 * abs(x)) ; }
key.foldChange = list(
  space = 'right',
  points = list(
    cex = key.foldChange.func(key.sizes),
    col = col.foldChange.func(key.sizes),
    pch = 19
  ),
  text = list(
    lab = as.character(key.sizes),
    cex = 1,
    adj = 1
  ),
  title = expression(bold(underline("Log"[2]*" FC"))),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)

# Difference in Medians
col.diffMed.func = function(i){
  cols = rep("darkorchid4", length(i))
  cols[i < 0] <- "turquoise3"
  cols[i == 0] <- "transparent"
  return(cols)
}

key.diffMed.sizes = seq(-10, 10, 5)
key.diffMed.func = function(x) { 0.1 + 0.3*abs(x) ; }
key.diffMed = list(
  space = 'right',
  points = list(
    cex = key.diffMed.func(key.diffMed.sizes),
    col = col.diffMed.func(key.diffMed.sizes),
    pch = 19
  ),
  text = list(
    lab = as.character(key.diffMed.sizes),
    cex = 1,
    adj = 1
  ),
  title = expression(bold(underline(Delta*" Medians"))),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)


# Hazard Ratio
col.hazardRatio.func = function(i){
  cols = rep("gold", length(i))
  cols[i < 0] <- "forestgreen"
  cols[i == 0] <- "transparent"
  return(cols)
}

key.hr.sizes = c(-2, -1, 0, 1, 2)
key.hazardRatio.func = function(x) { 0.1 + abs(x) ; }
key.hazardRatio = list(
  space = 'right',
  points = list(
    cex =  key.hazardRatio.func(key.hr.sizes),
    col = col.hazardRatio.func(key.hr.sizes),
    pch = 19
  ),
  text = list(
    lab = as.character(key.hr.sizes),
    cex = 1,
    adj = 1
  ),
  title = expression(bold(underline("Log"[2]*" HR"))),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

mut2num = structure(c(1, 2), names = c("Loss", "Gain"))
this.heatmap = function(i){
  plt.df = data.frame(
    "Mut" = mut2num[mutation.heatmap[[i]]],
    "Mut2" = mut2num[mutation.heatmap[[i]]]
  )
  plt = create.heatmap(
    plt.df,
    clustering.method = "none",
    same.as.matrix = T,
    # Axes
    ylab.label = firstup(names(mutation.heatmap)[i]),
    ylab.cex = 1.5,
    yaxis.lab = row.split[[i]],
    yaxis.cex = 1,
    yat = 1:nrow(plt.df),
    xaxis.tck = 0,
    yaxis.tck = 0,
    axes.lwd = 1,
    grid.row = TRUE,
    # Colours
    colour.scheme = c('dodgerblue', 'firebrick1'),
    at = c(0.5, 1.5, 2.5),
    print.colour.key = F
  )
  return(plt)
}

this.dotmap = function(i){
  # cat(i, "\n")

  # Function Assignment
  if(i %in% c(3, 6, 9)){
    key.func = key.hazardRatio.func
    col.func = col.hazardRatio.func
  } else if (i %in% c(2, 5, 8)){
    key.func = key.diffMed.func
    col.func = col.diffMed.func
  } else if (i %in% c(1, 4, 7)){
    key.func = key.foldChange.func
    col.func = col.foldChange.func
  }

  # Key Assignment
  if(i == 3){
    this.key = key.foldChange
  } else if (i == 6){
    this.key = key.diffMed
  } else if(i == 9){
    this.key = key.hazardRatio
  } else {
    this.key = NULL
  }

  plt = create.dotmap(
    #
    x = FC.list[[i]],
    bg.data = fdr.list[[i]],
    # Background colour
    colour.scheme = c('black', 'white'),
    at = c(0, 0.1, 1), # seq(0, 1, 0.1),
    bg.alpha = 1,
    # Spot functions
    spot.size.function = key.func,
    spot.colour.function = col.func,
    pch.border.col = 'white',
    # Scaling
    main.cex = 0,
    main = "",
    xaxis.lab = if(i %in% 7:9) col.split[[i-6]] else "",
    yaxis.lab = 0,
    xaxis.cex = if(i %in% 7:9) 1 else 0,
    # xaxis.rot = 90,
    yaxis.cex = 0,
    ylab.label = "",
    ylab.cex = 0,
    xlab.cex = 0,
    xaxis.tck = 0,
    yaxis.tck = 0,
    # lwd
    col.lwd = 1,
    row.lwd = 1,
    lwd = 1,
    # NA
    na.spot.size = 1,
    na.pch = 4,
    na.spot.size.colour = 'black',
    key = this.key,
    key.top = 1,
    colourkey = F, # if (i == 8) T else F,
    description = 'Dotmap created by BoutrosLab.plotting.general',
    resolution = 50
  )
  return(plt)
}

# Splitting
col.split = list(c("RNA", "Protein", "MeTPeak"), c("ISUP", "IDC", "T", "PSA", "Age", "Hypoxia"), c("BCR"))
df.compiled$Class = factor(df.compiled$Class, levels = c("writer", "reader", "eraser"))
row.split = split(df.compiled$Gene, df.compiled$Class)

split.panels = function(df){
  tmp = split(df, df.compiled$Class)
  tmp = lapply(tmp, function(i){
    return(list(
      i[,col.split[[1]], drop = F],
      i[,col.split[[2]], drop = F],
      i[,col.split[[3]], drop = F]))
  })
  unlist(tmp, recursive = F)
}

FC.list = split.panels(df.foldChange)
P.list = split.panels(df.Pvalue)
fdr.list = split.panels(df.Qvalue)

# Mutation Heatmap
mutation.heatmap = split(df.compiled$Mutation, df.compiled$Class)

# Making the plots
plts = lapply(1:12, function(i){
  cat(i, "\n")
  if(i %in% c(1, 5, 9)){
    this.heatmap(ceiling(i/4))
  } else if (i %in% 2:4) {
    this.dotmap(i - 1)
  } else if (i %in% 6:8) {
    this.dotmap( (i - 2))
  } else if (i %in% 10:12) {
    this.dotmap( (i - 3))
  }
})

filename = "~/figures/009_Summary.m6A.Enzyme.Associations.pdf"
pdf(filename, width = 11, height = 6)

create.multipanelplot(
  plot.objects = plts,
  layout.height = 3,
  layout.width = 4,
  plot.objects.heights = c(9, 18, 5),
  plot.objects.widths = c(4, 5.5, 10, 4.2),
  x.spacing = -4,
  y.spacing = -14
)

dev.off()
