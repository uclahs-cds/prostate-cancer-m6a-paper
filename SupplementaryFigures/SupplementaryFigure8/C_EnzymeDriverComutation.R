source("/cluster/home/helenzhu/code/snakemake_M6A_34_EnzymeMutations/000_HEADER.R")

library(BoutrosLab.plotting.general)
library(MetBrewer)

# Preamble ----------------------------------------------------------------
# This script runs a hypergeometric test for the co-mutation of m6A enzymes and drivers
# Gistic peaks are used to identify co-mutations which are expected
# Inspired by https://github.com/uclahs-cds/project-ProstateCancer-666pg/blob/640f2680e99e4a343e934b829226cc37e243d283/figure/driver_comparison/old_files/569pg_figure4a_driver_associations.R

# Loading Data ------------------------------------------------------------

# All Samples
all_samples = all_samples[!all_samples %in% exclude_samples]

# GISTIC Peaks
print(load("M_M6A_Enzyme_Mutations/Enzyme.GISTIC.Annotation.rsav"))
# [1] "gistic.results"  "non.gistic.muts" "muts.coords"  

# Creating Driver and Enzyme Matrices -------------------------------------

driver.matrix = mutations[,grep("reader|writer|eraser", colnames(mutations), invert = T)]

enzyme.matrix = mutations[,grep("reader|writer|eraser", colnames(mutations))]

# Analysis ----------------------------------------------------------------

# Hypergeometric test for co-mutation
comut.pval = matrix(
  0, 
  nrow = ncol(enzyme.matrix), 
  ncol = ncol(driver.matrix), 
  dimnames = list(
    colnames(enzyme.matrix), 
    colnames(driver.matrix))
  )

comut.effect.size = matrix(
  0, 
  nrow = ncol(enzyme.matrix), 
  ncol = ncol(driver.matrix), 
  dimnames = list(
    colnames(enzyme.matrix), 
    colnames(driver.matrix))
  )


for(i in rownames(comut.pval)){
  for(j in colnames(comut.pval)){
    
    combine.df = data.frame(
      "Enzyme" = enzyme.matrix[,i],
      "Driver" = driver.matrix[,j]
    )
    combine.df = combine.df[complete.cases(combine.df),]
    # This is testing for co-mutation (lower.tail = T for mutual exclusivity)
    comut.pval[i,j] <- phyper(
      q = sum(combine.df[,"Enzyme"] & combine.df[,"Driver"]), # number of white balls drawn
      m = sum(combine.df[,"Driver"]),# number of white balls
      n = nrow(combine.df)-sum(combine.df[,"Driver"]),# number of black balls
      k = sum(combine.df[,"Enzyme"]), # number of balls drawn
      lower.tail = F
    );
    
    observed = sum(combine.df[,"Enzyme"] & combine.df[,"Driver"])
    expected = sum(combine.df[,"Enzyme"])*sum(combine.df[,"Driver"])/nrow(combine.df)
    comut.effect.size[i, j] = observed - expected
  }
}

# Clustering --------------------------------------------------------------

row.order <- order.dendrogram(as.dendrogram(diana(comut.pval)))
comut.pval = comut.pval[row.order,]
col.order <- order.dendrogram(as.dendrogram(diana(t(comut.pval))))
comut.pval = comut.pval[,col.order]

comut.effect.size = comut.effect.size[row.order, col.order]

# For plotting
comut.pval = -log10(comut.pval)
comut.pval[comut.pval == Inf] <- -log10(.Machine$double.xmin)

# Doing a transpose
comut.pval = t(comut.pval)
comut.effect.size = t(comut.effect.size)

# Plotting ----------------------------------------------------------------

key.sizes = c(-30, -20, -10, 0, 10, 20, 30)
key.func = function(x) { 0.1 + log(abs(x)) ; }
key.effectSize = list(
  space = 'right',
  points = list(
    cex =  key.func(key.sizes),
    col = c("dodgerblue2", "dodgerblue2", "dodgerblue2", "transparent", "darkorange1", "darkorange1", "darkorange1"),
    pch = 19
  ),
  text = list(
    lab = as.character(key.sizes),
    cex = 0.8,
    adj = 1
  ),
  title = expression(bold(underline("Obs - Exp"))),
  cex.title = 0.8,
  padding.text = 5,
  background = 'white'
)


# Dotmap for association
association.dotmap = create.dotmap(
  #
  x = comut.effect.size,
  bg.data = comut.pval,
  # Background colour
  colour.scheme = c('white', 'black'),
  at = seq(0, 4, 1),
  bg.alpha = 1, 
  # Spot functions
  spot.size.function = key.func,
  # spot.colour.function = if(i %in% c(3, 6, 9)) col.hazardRatio.func else col.foldChange.func,
  pch.border.col = 'white',
  # Scaling
  main.cex = 0,
  main = "",
  xaxis.cex = 0,
  yaxis.cex = 0,
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
  key = key.effectSize,
  key.top = 1,
  colourkey = F,
  description = 'Dotmap created by BoutrosLab.plotting.general',
  resolution = 50
)


# Mutation Heatmaps -------------------------------------------------------

# Colours
standard.chrs = paste0("chr", c(1:22, "X", "Y"))
gistic.results$Gistic.Chr = factor(gistic.results$Gistic.Chr, levels = standard.chrs)
gistic.results = gistic.results[order(gistic.results$Gistic.Chr, gistic.results$Gistic.Start),]

encoding.col = c(
  met.brewer(name="Signac", n=length(unique(gistic.results$wide.peak.boundaries))),
  c("khaki3", "firebrick1", "dodgerblue")
)

encoding.num = structure(
  1:length(encoding.col),
  names = c(unique(gistic.results$wide.peak.boundaries), c("SSM", "Gain", "Loss"))
)

# Heatmap for GISTIC peak Enzymes
enzyme.covariate = colnames(comut.pval)
enzyme.covariate = do.call(rbind.data.frame, strsplit(enzyme.covariate, split = "[.]"))
colnames(enzyme.covariate) = c("Mutation", "Gene", "Class")
enzyme.covariate$GISTIC = gistic.results$wide.peak.boundaries[match(enzyme.covariate$Gene, gistic.results$Gene)]

# Heatmap for GISTIC peak Drivers
driver.covariate = rownames(comut.pval)
driver.covariate = do.call(rbind.data.frame, strsplit(driver.covariate, split = "[.]"))
colnames(driver.covariate) = c("Mutation", "Gene")
driver.covariate$Mutation[driver.covariate$Gene == "ETS"] <- "Loss"
driver.covariate$GISTIC = gistic.results$wide.peak.boundaries[match(driver.covariate$Gene, gistic.results$Gene)]
driver.covariate$GISTIC[driver.covariate$Mutation == "SSM"] <- NA
driver.covariate$GISTIC[driver.covariate$Gene == "ETS"] <- driver.covariate$GISTIC[driver.covariate$Gene == "ERG"]

# Encoding
enzyme.covariate$GISTIC = encoding.num[enzyme.covariate$GISTIC]
enzyme.covariate$Mutation = encoding.num[enzyme.covariate$Mutation]
driver.covariate$GISTIC = encoding.num[driver.covariate$GISTIC]
driver.covariate$Mutation = encoding.num[driver.covariate$Mutation]

# Making Plots ------------------------------------------------------------

enzyme.cov = create.heatmap(
  t(enzyme.covariate[,c("GISTIC", "Mutation")]),
  clustering.method = "none",
  same.as.matrix = T,
  # Axes
  xaxis.cex = 1,
  xat = 1:nrow(enzyme.covariate),
  xaxis.lab = enzyme.covariate$Gene,
  yaxis.cex = 0,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xlab.label = expression(bold("m"^"6"*"A Enzymes")),
  xlab.cex = 1.5,
  ylab.cex = 0,
  axes.lwd = 1,
  # Colour key
  colour.scheme = encoding.col,
  at = seq(0.5, length(encoding.col) + 0.5, 1),
  total.colours = length(encoding.col) + 1,
  fill.colour = 'white', 
  # colourkey.cex = 1,
  # colourkey.labels.at = 1:length(encoding.col),
  # colourkey.labels = names(encoding.num),
  print.colour.key = F
)

driver.cov = create.heatmap(
  driver.covariate[,c("Mutation", "GISTIC")],
  clustering.method = "none",
  same.as.matrix = T,
  # Axes
  yaxis.cex = 1,
  yat = 1:nrow(driver.covariate),
  yaxis.lab = driver.covariate$Gene,
  xaxis.cex = 0,
  xaxis.tck = 0,
  yaxis.tck = 0,
  ylab.label = "Prostate Cancer Drivers",
  xlab.cex = 0,
  ylab.cex = 1.5,
  axes.lwd = 1,
  # Colour key
  colour.scheme = encoding.col,
  at = seq(0.5, length(encoding.col) + 0.5, 1),
  total.colours = length(encoding.col) + 1,
  fill.colour = 'white', 
  # colourkey.cex = 1,
  # colourkey.labels.at = 1:length(encoding.col),
  # colourkey.labels = names(encoding.num),
  print.colour.key = F
)


# Legend
colfunc <- colorRampPalette(c("white", "black"))
covariate.legend <- list(
  legend = list(
    colours =  met.brewer(name="Signac", n=length(unique(gistic.results$wide.peak.boundaries))),
    labels = unique(gistic.results$wide.peak.boundaries),
    title = expression(bold(underline('GISTIC Peaks'))),
    lwd = 0.5
  ),
  legend = list(
    colours = c("khaki3", "firebrick1", "dodgerblue"),
    labels =  c("SSM", "Gain", "Loss"),
    title = expression(bold(underline('Mutations'))),
    lwd = 0.5
  ),
  legend = list(
    colours = colfunc(4),
    labels = c(
      expression("P < 10"^"-1"),
      expression("P < 10"^"-2"),
      expression("P < 10"^"-3"),
      expression("P < 10"^"-4")
    ),
    title = expression(bold(underline('P value'))),
    lwd = 0.5
  )
)

side.legend <-legend.grob(
  legends = rev(covariate.legend),
  label.cex = 0.8,
  title.cex = 0.8,
  title.just = 'left',
  title.fontface = 'bold',
  size = 2
)

filename = "~/figures/020_EnzymeDriverComutationPlot.pdf"
pdf(filename, width = 9, height = 9)

create.multipanelplot(
  plot.objects = list(driver.cov, association.dotmap, enzyme.cov),
  layout.height = 2,
  layout.width = 2,
  plot.objects.heights = c(5, 1.8),
  plot.objects.widths = c(2, 5),
  y.spacing = -2.1,
  x.spacing = -2.5,
  layout.skip = c(F, F,
                  T, F),
  legend = list(
    right = list(
      x = 0.8,
      y = 1,
      fun = side.legend
    )
  ),
  right.padding = 1
)
dev.off()