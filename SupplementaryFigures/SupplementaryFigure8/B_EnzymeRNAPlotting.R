source("/cluster/home/helenzhu/code/snakemake_M6A_34_EnzymeMutations/000_HEADER.R")


# Preamble ----------------------------------------------------------------
# This makes a dotmap of the mutation associations
# RNA
# Protein
# Number of Peaks
# Clinical Variables

# Loading Data ------------------------------------------------------------

# RNA Results
print(load("M_M6A_Enzyme_Mutations/RNA.m6a.enzymes.clinical.data.rsav"))
# [1] "rna"         "rna.results"

# Protein Results
print(load("M_M6A_Enzyme_Mutations/Protein.m6a.enzymes.clinical.data.rsav"))
# [1] "protein"         "protein.results"

# Creating a Dotmap -------------------------------------------------------

col.foldChange.func = function(i){
  cols = rep("darkorange1", length(i))
  cols[i < 0] <- "dodgerblue2"
  cols[i == 0] <- "transparent"
  return(cols)
}

key.sizes = c(-0.5, -0.1, 0, 0.1, 0.5)
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
  title = expression(bold(underline(rho))),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)

col.hazardRatio.func = function(i){
  cols = rep("gold", length(i))
  cols[i < 0] <- "forestgreen"
  cols[i == 0] <- "transparent"
  return(cols)
}

key.hr.sizes = c(-1.5, -1, 0, 1, 1.5)
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

enzyme.levs = c("writer", "reader", "eraser")
this.dotmap = function(i, FC.list, fdr.list, row.split){
  # cat(i, "\n")
  plt = create.dotmap(
    #
    x = FC.list[[i]],
    bg.data = fdr.list[[i]],
    # Background colour
    colour.scheme = c('black', 'white'),
    at = c(0, 0.1, 1), # seq(0, 1, 0.1),
    bg.alpha = 1,
    # Spot functions
    spot.size.function = if(i %in% c(3, 6, 9)) key.hazardRatio.func else key.foldChange.func,
    spot.colour.function = if(i %in% c(3, 6, 9)) col.hazardRatio.func else col.foldChange.func,
    pch.border.col = 'white',
    # Scaling
    main.cex = 0,
    main = "",
    xaxis.lab = if(i %in% 7:9) col.split[[i-6]] else "",
    xaxis.cex = if(i %in% 7:9) 1 else 0,
    # xaxis.rot = 90,
    yaxis.cex = if(i %in% c(1, 4, 7)) 1 else 0,
    yaxis.lab = if(i %in% c(1, 4, 7)) row.split[[ceiling(i/3)]] else "",
    ylab.label = if(i %in% c(1, 4, 7)) firstup(enzyme.levs[ceiling(i/3)]) else "",
    ylab.cex = if(i %in% c(1, 4, 7)) 1 else 0,
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
    key = if(i == 3) key.foldChange else if (i == 6) key.hazardRatio else NULL,
    key.top = 1,
    colourkey = F, # if (i == 8) T else F,
    description = 'Dotmap created by BoutrosLab.plotting.general',
    resolution = 50
  )
  return(plt)
}

# Splitting
col.split = list(c("MeTPeak"), c("ISUP", "IDC", "T", "PSA", "Age", "Hypoxia"), c("BCR"))
split.panels = function(df, enzyme.class){
  tmp = split(df, enzyme.class)
  tmp = lapply(tmp, function(i){
    return(list(
      i[,col.split[[1]], drop = F],
      i[,col.split[[2]], drop = F],
      i[,col.split[[3]], drop = F]))
  })
  unlist(tmp, recursive = F)
}



# Making Plots ------------------------------------------------------------

generate.plots = function(df.compiled){

  # df.compiled = rna.results

  # Q value
  df.Qvalue = df.compiled[,grep(".fdr", colnames(df.compiled), value = T, fixed = T)]
  colnames(df.Qvalue) = gsub(".fdr", "", colnames(df.Qvalue), fixed = T)

  # Fold Change
  df.foldChange = df.compiled[,grep(".effectSize", colnames(df.compiled), value = T, fixed = T)]
  df.foldChange$BCR.effectSize = log2(df.foldChange$BCR.effectSize)
  colnames(df.foldChange) = gsub(".effectSize", "", colnames(df.foldChange), fixed = T)

  missing.cols = setdiff(colnames(df.Qvalue), colnames(df.foldChange))
  df.fc.missing = matrix(0, ncol = length(missing.cols), nrow = nrow(df.foldChange), dimnames = list(rownames(df.foldChange), missing.cols))
  df.foldChange = cbind(df.foldChange, df.fc.missing)
  df.foldChange = df.foldChange[,colnames(df.Qvalue)]

  # Check
  all(colnames(df.foldChange) == colnames(df.Qvalue))

  # Row names
  rnames = df.compiled$Gene

  # Col names
  cnames = colnames(df.Qvalue)

  # Enzyme class
  enzyme.class = enzymes$Class[match(rnames, enzymes$GeneName)]
  enzyme.class = factor(enzyme.class, levels = c("writer", "reader", "eraser"))
  row.split = split(rnames, enzyme.class)

  # Splitting Panels
  FC.list = split.panels(df.foldChange, enzyme.class)
  fdr.list = split.panels(df.Qvalue, enzyme.class)

  # Making the plots
  plts = lapply(1:9, function(i){
    cat(i, "\n")
    this.dotmap(i, FC.list, fdr.list, row.split)
  })

  plts
}

# RNA plots
rna.plts = generate.plots(rna.results)
filename = "~/figures/033_RNA.m6A.Enzyme.Associations.pdf"
pdf(filename, width = 8, height = 6)

create.multipanelplot(
  plot.objects = rna.plts,
  layout.height = 3,
  layout.width = 3,
  plot.objects.heights = c(14, 20, 7),
  plot.objects.widths = c(4, 7.5, 3.4),
  x.spacing = -3.5,
  y.spacing = -15
)

dev.off()

# Protein plots
protein.plts = generate.plots(protein.results)
filename = "~/figures/033_Protein.m6A.Enzyme.Associations.pdf"
pdf(filename, width = 8, height = 6)

create.multipanelplot(
  plot.objects = protein.plts,
  layout.height = 3,
  layout.width = 3,
  plot.objects.heights = c(14, 20, 7),
  plot.objects.widths = c(4, 7.5, 3.4),
  x.spacing = -3.5,
  y.spacing = -15
)

dev.off()
