source("/cluster/home/helenzhu/code/snakemake_M6A_8_PipelineRebuild/Rscripts/000_HEADER.R")

library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(valr)

# Preamble ----------------------------------------------------------------
# This makes plots for genes
# TODO:
# 1. Fix scale for plots
# 2. Add an inset for AR
# 3. Update plots on GDrive

# AR research
# APPRIS database on primary transcript: https://appris.bioinfo.cnio.es/#/database/id/homo_sapiens/ENSG00000169083?as=hg38&sc=ensembl&ds=e103v45
# Nature reviews cancer on resistance: https://www.nature.com/articles/nrc4016

# Loading Data ------------------------------------------------------------

# Filtering Samples
all_samples = all_samples[!all_samples %in% exclude_samples]

# MALAT1, PTEN, MYC, FOXA1, TP53, AR, SPOP, RB1, NKX3-1
genes = c("MALAT1", "PTEN", "MYC", "FOXA1", "TP53", "AR", "SPOP", "RB1", "NKX3-1")
genes = gene_info[gene_info$Name %in% genes,]
genes = genes[genes$GeneID != "ENSG00000284792.1",]

# AR principal isoform
ar_principal <- "ENST00000374690.9"
ar_coords <- c(67723839, 67730619) # Last UTR of ENSG00000169083

# This is the last four exons/cds that make up the ligand binding domain, which is frequently mutated in prostate cancer
# c(67717478, 67723838) 

# Loading gene annotations
gtf.file = "GRCh38.p13_Build34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf"
gtf = read_gtf(gtf.file, zero_based = F)

# Creating a GRanges Object
genes.gr = gtf[gtf$gene_id %in% genes$GeneID &
                 gtf$type == "gene",
               c("chrom", "start", "end", "strand", "gene_id", "gene_name")]
genes.gr = makeGRangesFromDataFrame(genes.gr, keep.extra.columns = T)

# Gene Coordinates
# & gtf$type %in% c("exon", "UTR")
genes.coords = gtf[gtf$gene_id %in% genes$GeneID,
                   c("chrom", "start", "end", "strand", "gene_id", "gene_name", "type")]
genes.coords = data.frame(genes.coords)

# Calculate coverage
calculate.coverage = function(f.cov, r.cov, genes.gr){
  histogram.coverage = vector("list", length(genes.gr))
  names(histogram.coverage) = genes.gr$gene_id
  for(x in seq_along(genes.gr)){
    gr = genes.gr[x]
    bins = unlist(GenomicRanges::tile(x = gr, width = 1))
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(f.cov)
    if(as.character(strand(gr)) == "+"){
      cvg = binnedAverage(
        bins = bins,
        numvar = f.cov,
        varname = "cvg")
    } else {
      cvg = binnedAverage(
        bins = bins,
        numvar = r.cov,
        varname = "cvg")
    }
    histogram.coverage[[gr$gene_id]] <- cvg$cvg
  }
  return(histogram.coverage)
}

# IP Bigwigs
# ip.coverage = lapply(all_samples, function(sample.id){
#   cat(sample.id, "\n")
#   # Forward BigWig
#   f.bw.file = paste0("BigWigs/", sample.id, "_IP.forward.bigWig")
#   f.bw = valr::read_bigwig(f.bw.file, set_strand = "+")
#   f.bw$start = f.bw$start + 1
#   f.gr = makeGRangesFromDataFrame(f.bw, keep.extra.columns = T)
#   f.gr = f.gr[f.gr$score > 0]
#   f.cov = coverage(f.gr, weight = "score")
# 
#   # Reverse BigWig
#   r.bw.file = paste0("BigWigs/", sample.id, "_IP.reverse.bigWig")
#   r.bw = valr::read_bigwig(r.bw.file, set_strand = "-")
#   r.bw$start = r.bw$start + 1
#   r.gr = makeGRangesFromDataFrame(r.bw, keep.extra.columns = T)
#   r.gr = r.gr[r.gr$score > 0]
#   r.cov = coverage(r.gr, weight = "score")
# 
#   # Calculate coverage
#   calculate.coverage(f.cov, r.cov, genes.gr)
# })

# Input Bigwigs
# input.coverage = lapply(all_samples, function(sample.id){
#   cat(sample.id, "\n")
#   # Forward BigWig
#   f.bw.file = paste0("BigWigs/", sample.id, "_Input.forward.bigWig")
#   f.bw = valr::read_bigwig(f.bw.file, set_strand = "+")
#   f.bw$start = f.bw$start + 1
#   f.gr = makeGRangesFromDataFrame(f.bw, keep.extra.columns = T)
#   f.gr = f.gr[f.gr$score > 0]
#   f.cov = coverage(f.gr, weight = "score")
# 
#   # Reverse BigWig
#   r.bw.file = paste0("BigWigs/", sample.id, "_Input.reverse.bigWig")
#   r.bw = valr::read_bigwig(r.bw.file, set_strand = "-")
#   r.bw$start = r.bw$start + 1
#   r.gr = makeGRangesFromDataFrame(r.bw, keep.extra.columns = T)
#   r.gr = r.gr[r.gr$score > 0]
#   r.cov = coverage(r.gr, weight = "score")
# 
#   # Calculate coverage
#   calculate.coverage(f.cov, r.cov, genes.gr)
# })

# names(ip.coverage) = all_samples
# names(input.coverage) = all_samples
# save(ip.coverage, input.coverage, file = "prostate.gene.coverage.rsav")

# Loading saved data
print(load("prostate.gene.coverage.rsav"))
# [1] "ip.coverage"    "input.coverage"


# Extracting Genes --------------------------------------------------------

ip.genes = lapply(genes$GeneID, function(i) {
  do.call(rbind, lapply(ip.coverage, '[[', i))
})
names(ip.genes) = genes$Name

input.genes = lapply(genes$GeneID, function(i) {
  do.call(rbind, lapply(input.coverage, '[[', i))
})
names(input.genes) = genes$Name

# Plotting ----------------------------------------------------------------

# Coordinates refer to a subset of the gene - a vector of length 2
# c(start_coordinate, end_coordinate)
plot_gene <- function(gene_id, coordinates = NULL, max.data = 4000, max.data.by = 1000, inset.lines = NULL, filename){

  # gene_id = "MYC"
  cat(gene_id, "\n")

  # Plotting Data
  if(!is.null(coordinates)){
    gene.start = genes.coords[
      genes.coords$gene_name == gene_id & genes.coords$type == "gene", 
      "start"
    ]
    coordinates <- coordinates - gene.start + 1
    ip.data = ip.genes[[gene_id]][,coordinates[1]:coordinates[2]]
    input.data = input.genes[[gene_id]][,coordinates[1]:coordinates[2]]
  } else {
    ip.data = ip.genes[[gene_id]]
    input.data = input.genes[[gene_id]]
  }

  # Plotting parameters
  rect.range = -0.1*max.data
  
  # Polygonplot plotting data
  pg.data = data.frame(
    "x" = rep(1:ncol(ip.data), 2),
    "min" = c(apply(ip.data, 2, min), apply(input.data, 2, min)),
    "max" = c(apply(ip.data, 2, max), apply(input.data, 2, max)),
    "median" = c(apply(ip.data, 2, median), apply(input.data, 2, median)),
    "groups" = c(rep("IP", ncol(ip.data)), rep("Input", ncol(ip.data)))
  )

  # Gene Annotations
  rect.data = genes.coords[genes.coords$gene_name == gene_id,]
  gene.start = rect.data$start[rect.data$type == "gene"]
  rect.data$start = rect.data$start - gene.start + 1
  rect.data$end = rect.data$end - gene.start + 1
  rect.data = rect.data[rect.data$type %in% c("UTR", "exon"), c('start', 'end', 'type')]
  rect.data$bottom = ifelse(rect.data$type == "UTR", 0.7*rect.range, 0.9*rect.range)
  rect.data$top = ifelse(rect.data$type == "UTR", 0.3*rect.range, 0.1*rect.range)
  rect.data = rbind(rect.data, data.frame('start' = 1, 'end' = ncol(ip.data), 'type' = 'gene', 'bottom' = 0.55*rect.range, 'top' = 0.45*rect.range))
  rect.data$col <- "darkblue"
  
  if(!is.null(inset.lines)){
    gene.start = genes.coords[
      genes.coords$gene_name == gene_id & genes.coords$type == "gene", 
      "start"
    ]
    inset.lines <- inset.lines - gene.start + 1
    rect.data <- rbind(
      rect.data, 
      data.frame('start' = inset.lines[1], 'end' = inset.lines[2], type = 'inset', 'bottom' = 0, 'top' = max.data, 'col' = 'pink')
    )
  }
  
  if(!is.null(coordinates)){
    # Find overlapping annotations
    rect.data <- rect.data[
      (rect.data$start >= coordinates[1] & rect.data$start <= coordinates[2]) | # Start overlapps coordinates
        (rect.data$end >= coordinates[1] & rect.data$end <= coordinates[2]) | 
        (rect.data$start <= coordinates[1] & rect.data$end >= coordinates[2]),
      ]
    # Shift coordinates of annotations
    rect.data$start <- rect.data$start - coordinates[1] + 1
    rect.data$end <- rect.data$end - coordinates[1] + 1
    rect.data$start <- ifelse(rect.data$start <= 1, 1, rect.data$start)
    rect.data$end <- ifelse(rect.data$end >= ncol(ip.data), ncol(ip.data), rect.data$end)
  }
  
  # Formatting plot
  if (max.data <= 15000){
    max.data.scale = 10^3
    ylab_label = expression(bold("Coverage (10"^"3"*")"))
  } else if (max.data > 15000 && max.data <= 150000){
    max.data.scale = 10^4
    ylab_label = expression(bold("Coverage (10"^"4"*")"))
  } else if (max.data > 150000){
    max.data.scale = 10^5
    ylab_label = expression(bold("Coverage (10"^"5"*")"))
  }
  
  this.yat = seq(0, max.data, max.data.by) # Changed max.data to 4000
  this.yaxis.lab = this.yat/max.data.scale

  # create the simple plot
  pg.plt = create.polygonplot(
    formula = NA ~ x,
    data = pg.data,
    max = pg.data$min,
    min = pg.data$max,
    groups = pg.data$groups,
    # Main
    # main = gene_id,
    # main.just = 'left',
    # main.x = 0.11,
    # main.y = -0.4,
    # main.cex = 1.5,
    # Adding lines for inset
    # abline.v = inset.lines,
    # abline.col = 'red',
    # Adding gene name as text
    add.text = if(!is.null(coordinates)) FALSE else TRUE,
    text.labels = gene_id, # pad_nchar(gene_id),
    text.x = 0.5*ncol(ip.data),
    text.y = 0.875*max.data,
    text.col = 'black',
    text.cex = 4,
    text.fontface = 'bold',
    # adjust axes limits
    xlimits = c(1, ncol(ip.data)),
    ylimits = c(rect.range, max.data),
    ylab.cex = if(!is.null(coordinates)) 4 else 2,
    xlab.cex = 0,
    yaxis.cex = if(!is.null(coordinates)) 4 else 2,
    xaxis.cex = 0,
    yaxis.tck = 0,
    xaxis.tck = 0,
    yat = this.yat,
    yaxis.lab = this.yaxis.lab,
    ylab.label = ylab_label,
    # add fill colour
    col = c("darkorchid4", "seagreen3"),
    border.col = c("darkorchid4", "seagreen3"),
    alpha = 0.2,
    # add middle line
    add.median = TRUE,
    median = pg.data$median,
    median.lty = 1,
    median.lwd = 1,
    median.col = c("darkorchid4", "seagreen3"),
    # Add rectangles
    add.rectangle = T,
    xleft.rectangle = rect.data$start,
    xright.rectangle = rect.data$end,
    ybottom.rectangle = rect.data$bottom,
    ytop.rectangle = rect.data$top,
    col.rectangle = rect.data$col, # 'darkblue',
    alpha.rectangle = 1
  )

  log.ip.data = log1p(ip.data)
  log.input.data = log1p(input.data)
  max.log.data = max(rbind(log.input.data, log.ip.data))+0.1

  # Creating a heatmap for IP
  ip.hm = create.heatmap(
    log.ip.data,
    # Clustering
    clustering.method = "ward.D2",
    cluster.dimensions = "row",
    rows.distance.method = "euclidean", # correlation",
    plot.dendrograms = 'none',
    # Plotting
    same.as.matrix = T,
    xaxis.cex = 0,
    yaxis.cex = 0,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xat = seq(1, ncol(ip.data), length.out = 4),
    yaxis.lab = rownames(ip.data),
    xaxis.lab = seq(1, ncol(ip.data), length.out = 4),
    axes.lwd = 1,
    ylab.label = if(is.null(coordinates)) "IP" else "",
    ylab.cex = 2,
    # Colours
    colour.scheme = c('white', 'seagreen3'),
    at = seq(0, max.log.data, 0.1),
    print.colour.key = F
  )

  # Creating a heatmap for Input
  input.hm = create.heatmap(
    log.input.data,
    # Clustering
    clustering.method = "ward.D2",
    cluster.dimensions = "row",
    rows.distance.method = "euclidean", # "correlation",
    plot.dendrograms = 'none',
    # Plotting
    same.as.matrix = T,
    xaxis.cex = 0,
    yaxis.cex = 0,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xat = seq(1, ncol(ip.data), length.out = 4),
    yaxis.lab = rownames(ip.data),
    xaxis.lab = seq(1, ncol(ip.data), length.out = 4),
    axes.lwd = 1,
    ylab.label = if(is.null(coordinates)) "Input" else "",
    ylab.cex = 2,
    # Colours
    colour.scheme = c('white', 'darkorchid4'),
    at = seq(0, max.log.data, 0.1),
    print.colour.key = F
  )

  # Creating a multipanelplot
  mplt = create.multipanelplot(
    plot.objects = list(pg.plt, ip.hm, input.hm),
    plot.objects.heights = c(3, 1, 1),
    y.spacing = -2,
    ylab.axis.padding = -0.5,
    # Printing plot
    filename = filename,
    height = 9,
    width = 8,
    resolution = 1000,
    size.units = 'in'
  )

  return(mplt)
}

# Applying plot_gene to all genes

# MYC
filename = "~/figures/400_GeneProfilePlots.MYC.png"
plot_gene("MYC", filename = filename)

# TP53
filename = "~/figures/400_GeneProfilePlots.TP53.png"
plot_gene("TP53", filename = filename)

# PTEN
filename = "~/figures/400_GeneProfilePlots.PTEN.png"
plot_gene("PTEN", filename = filename)

# AR
filename = "~/figures/400_GeneProfilePlots.AR.png"
plot_gene("AR", inset.lines = ar_coords, filename = filename)

# AR Inset
filename = "~/figures/400_GeneProfilePlots.AR.Inset.png"
plot_gene("AR", coordinates = ar_coords, max.data = 4000, max.data.by = 4000, filename = filename)

# MALAT1
filename = "~/figures/400_GeneProfilePlots.MALAT1.png"
plot_gene("MALAT1", max.data = 1200000, max.data.by = 400000, filename = filename)

# FOXA1
filename = "~/figures/400_GeneProfilePlots.FOXA1.png"
plot_gene("FOXA1", max.data = 150000, max.data.by = 50000, filename = filename)

# RB1
filename = "~/figures/400_GeneProfilePlots.RB1.png"
plot_gene("RB1", max.data = 1000, max.data.by = 500, filename = filename)

# NKX3-1
filename = "~/figures/400_GeneProfilePlots.NKX3-1.png"
plot_gene("NKX3-1", max.data = 45000, max.data.by = 20000, filename = filename)