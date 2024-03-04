source("/cluster/home/helenzhu/code/snakemake_M6A_9_PeakMotifs/Rscripts/000_HEADER.R")

# library(MeRIPtools)
# library(RADAR)

# Preamble ----------------------------------------------------------------
# Comparison of Peak Callers
# Comparison to the REPIC database

# Loading Data ------------------------------------------------------------

# Samples
all_samples = all_samples[!all_samples %in% exclude_samples]

# MeTPeak
metpeak_peaks = do.call(rbind, lapply(all_samples, function(this_sample){
  peaks = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/F_MeTPeak_Peaks_SS/", this_sample, "/peak.bed6"),
                     header = F, sep = "\t", stringsAsFactors = F)
  peaks$sample = this_sample
  peaks
}))

# ExomePeak
exomepeak_peaks = do.call(rbind, lapply(all_samples, function(this_sample){
  peaks = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/G_exomePeak_Peaks_SS/", this_sample, "/peak.bed6"),
                     header = F, sep = "\t", stringsAsFactors = F)
  peaks$sample = this_sample
  peaks
}))

# Number of peaks
metpeak_num = do.call(rbind, lapply(all_samples, function(this_sample){
  peaks = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/F_MeTPeak_Peaks_SS/", this_sample, "/peak.bed"),
                     header = F, sep = "\t", stringsAsFactors = F)
  data.frame("sample" = this_sample, "MeTPeak_peaks" = nrow(peaks), stringsAsFactors = F)
}))

exomepeak_num = do.call(rbind, lapply(all_samples, function(this_sample){
  peaks = read.table(paste0("/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/G_exomePeak_Peaks_SS/", this_sample, "/peak.bed"),
                     header = F, sep = "\t", stringsAsFactors = F)
  data.frame("sample" = this_sample, "exomePeak_peaks" = nrow(peaks), stringsAsFactors = F)
}))

# REPIC humans
repic = read.table("/cluster/home/helenzhu/Cluster_Helen/REPIC/m6a.sites.species.human.hg38.bed",
                   header = F, sep = "\t", stringsAsFactors = F)[1:5]

# Sequencing Stats
print(load("Summary_Tables/hg38.STAR.output.rsav"))

# Peak Similarity ---------------------------------------------------------

# JACCARD
jaccard_peaks = function(peaks1, peaks2){

  # Filtering Chrs
  chrs = paste0("chr", c(1:22, "X", "Y"))
  peaks1 = peaks1[peaks1[,1] %in% chrs,]
  peaks2 = peaks2[peaks2[,1] %in% chrs,]

  peaks1_gr = GRanges(
    seqnames = peaks1[,1],
    IRanges(start = peaks1[,2],
            end = peaks1[,3])
  )
  peaks2_gr = GRanges(
    seqnames = peaks2[,1],
    IRanges(start = peaks2[,2],
            end = peaks2[,3])
  )

 jaccard_index = sum(width(intersect(peaks1_gr, peaks2_gr)))/sum(width(union(peaks1_gr, peaks2_gr)))
 return(jaccard_index)

}

# Simpson
simpson_peaks = function(peaks1, peaks2){

  # Filtering Chrs
  chrs = paste0("chr", c(1:22, "X", "Y"))
  peaks1 = peaks1[peaks1[,1] %in% chrs,]
  peaks2 = peaks2[peaks2[,1] %in% chrs,]

  peaks1_gr = GRanges(
    seqnames = peaks1[,1],
    IRanges(start = peaks1[,2],
            end = peaks1[,3])
  )
  peaks2_gr = GRanges(
    seqnames = peaks2[,1],
    IRanges(start = peaks2[,2],
            end = peaks2[,3])
  )

  nm = sum(width(intersect(peaks1_gr, peaks2_gr)))
  dm = min(c(sum(width(reduce(peaks1_gr))), sum(width(reduce(peaks2_gr)))))
  simpson_index = nm/dm

  return(simpson_index)
}

peak_stats = do.call(rbind, lapply(all_samples, function(this_sample){

  cat(this_sample, "\n")
  metpeak = metpeak_peaks[metpeak_peaks$sample == this_sample,]
  exomepeak = exomepeak_peaks[exomepeak_peaks$sample == this_sample,]

  data.frame(
    "Sample" = this_sample,
    "Jaccard_MeTPeakExomePeak" = jaccard_peaks(metpeak, exomepeak),
    "Simpson_MeTPeakExomePeak" = simpson_peaks(metpeak, exomepeak),
    stringsAsFactors = F
  )

}))

# Jaccard-Simpson plot
make_scatter = function(col1, col2, this_main){

  scatter_df = data.frame(
    "x" = peak_stats[,col1]*100,
    "y" = peak_stats[,col2]*100,
    stringsAsFactors = F
  )

  tmp = create.scatterplot(
    y ~ x,
    data = scatter_df,
    # main = this_main,
    # main.cex = 1,
    xlimits = c(0, 100),
    ylimits = c(0, 100),
    xat = seq(0, 100, 25),
    yat = seq(0, 100, 25),
    xaxis.tck = 0,
    yaxis.tck = 0,
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    ylab.label = "Simpson Index (%)",
    xlab.label = "Jaccard Index (%)",
    ylab.cex = 1.5,
    xlab.cex = 1.5,
    alpha = 0.2,
    cex = 1.5
  )

  return(tmp)
}

filename = "figures/100_JaccardSimpsonPlot.pdf"
pdf(filename, width = 4, height = 3.4)

make_scatter(
  col1 = "Jaccard_MeTPeakExomePeak",
  col2 = "Simpson_MeTPeakExomePeak",
  this_main = "XX")

dev.off()
system(paste0("cp ", filename, " ~/figures"))

# todo = c("MeTPeakExomePeak", "MeTPeakMACS2", "MeTPeakBinomial", "ExomePeakMACS2", "ExomePeakBinomial", "MACS2Binomial", "AllSingleBinomial")
# this_main = c("MeTPeak & ExomePeak", "MeTPeak & MACS2", "MeTPeak & Binomial", "ExomePeak & MACS2", "ExomePeak & Binomial", "MACS2 & Binomial", "MACS2,ExomePeak,MeTPeak & Binomial")
# todo = c("MeTPeakExomePeak")
# this_main = c("MeTPeak & ExomePeak")
# todo = data.frame("col1" = paste0("Jaccard_", todo),
#                   "col2" = paste0("Simpson_", todo),
#                   "main" = this_main,
#                   stringsAsFactors = F)
# this_plots = lapply(1:nrow(todo), function(i) make_scatter(col1 = todo$col1[i], col2 = todo$col2[i], this_main = todo$main[i]))

# create.multipanelplot(
#   plot.objects = this_plots,
#   layout.height = 1, # 3,
#   layout.width = 1,
#   y.spacing = 1,
#   x.spacing = 1,
#   ylab.label = "Simpson Index",
#   xlab.label = "Jaccard Index",
#   ylab.cex = 0,
#   xlab.cex = 0
# )

# Robustness to Read Depth ------------------------------------------------

ip_num = star_output[,c("Sample", "UniquelyMapped.Count.IP")]
colnames(ip_num) = c("sample", "IP_reads")
input_num = star_output[,c("Sample", "UniquelyMapped.Count.Input")]
colnames(input_num) = c("sample", "Input_reads")

peaknum_reads = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sample", all.x = TRUE),
                       list(metpeak_num, exomepeak_num, ip_num, input_num))

make_reads_scatter = function(col1, col2, this_ylab, this_xlab, this_yaxis.cex, this_xaxis.cex){

  scatter_df = data.frame(
    "x" = peaknum_reads[,col1]/10^6,
    "y" = peaknum_reads[,col2]/1000,
    stringsAsFactors = F
  )

  tmp = create.scatterplot(
    y ~ x,
    data = scatter_df,
    xaxis.cex = this_xaxis.cex,
    yaxis.cex = this_yaxis.cex,
    ylab.label = this_ylab,
    xlab.label = this_xlab,
    xlab.cex = 1.5,
    ylab.cex = 1.5,
    # xaxis.rot = 90,
    ylimits = c(0, 30),
    xlimits = c(0, 150),
    yat = seq(0, 25, 10),
    xat = seq(0, 140, 50),
    xaxis.tck = 0,
    yaxis.tck = 0,
    # Point features
    col = 'black',
    alpha = 0.2,
    cex = 1.5,
    # pch = 1,
    # alpha = 0.5,
    legend = list(
      inside = list(
        fun = draw.key,
        args = list(
          key = get.corr.key(
            x = scatter_df$x,
            y = scatter_df$y,
            label.items = c('spearman','spearman.p'),
            alpha.background = 0,
            key.cex = 1.5
          )
        ),
        x = 0,
        y = 1,
        corner = c(0,1)
      )
    )
  )

  return(tmp)
}

todo = expand.grid(c("MeTPeak_peaks", "exomePeak_peaks"),
                   c("IP_reads", "Input_reads"),
                   stringsAsFactors = F)
todo$xlab = c("", "", "MeTPeak", "exomePeak")
todo$ylab = c("IP", "", "Input", "")
todo$xaxis.cex = c(0, 0, 1.5, 1.5)
todo$yaxis.cex = c(1.5, 0, 1.5, 0)

this_plots = lapply(1:nrow(todo), function(i)
  make_reads_scatter(
    col1 = todo$Var2[i],
    col2 = todo$Var1[i],
    this_ylab = todo$ylab[i],
    this_xlab = todo$xlab[i],
    this_yaxis.cex = todo$yaxis.cex[i],
    this_xaxis.cex = todo$xaxis.cex[i]
    )
  )

filename = "figures/100_Peaks.vs.LibrarySize.pdf"
pdf(filename, width = 7.5, height = 7.5)

create.multipanelplot(
  plot.objects = this_plots,
  layout.height = 2,
  layout.width = 2,
  plot.objects.heights = c(1, 1.1),
  plot.objects.widths = c(1.15, 1),
  y.spacing = -1,
  x.spacing = 0,
  ylab.label = "# Peaks (per thousand)",
  xlab.label = "# Uniquely mapped reads (per million)",
  ylab.cex = 1.5,
  xlab.cex = 1.5
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))

# bold(paste("# Peaks (10"["3"]*")"))
# bold(paste("# Uniquely mapped reads (10"["6"]*")"))

# Paper Figures
correlation_results = do.call(rbind, lapply(1:nrow(todo), function(i){

  tmp = cor.test(
    peaknum_reads[,todo$Var1[i]],
    peaknum_reads[,todo$Var2[i]],
    method = "spearman")
  pval = tmp$p.value
  rho = tmp$estimate
  data.frame(
    "Peaks" = todo$Var1[i],
    "Library" = todo$Var2[i],
    "P" = pval,
    "Rho" = rho,
    stringsAsFactors = F
  )
}))

correlation_results$Peaks = gsub("_peaks", "", correlation_results$Peaks)
correlation_results$Library = gsub("_reads", "", correlation_results$Library)

pval = acast(data = correlation_results, Peaks ~ Library, value.var = "P")
pval = -log10(pval)
rho = acast(data = correlation_results, Peaks ~ Library, value.var = "Rho")


changeSciNot <- function(n) {
  output <- formatC(n, digits = 2, format = "e")
  num = strsplit(output, split = "e")[[1]][1]
  exponent = strtoi(gsub("-0", "-", strsplit(output, split = "e")[[1]][2]))
  if(num == "1.00"){
    output = bquote("" ~ 10^.(exponent))
  }else{
    output = bquote("" ~ .(num) ~ x ~ 10^.(exponent))
  }
  output
}


key.sizes = c(0.2, 0.3, 0.4, 0.5)
spot.size.function = function(x) 0.1 + (8 * abs(x))

filename = "figures/100_Peaks.vs.LibrarySize.Paper.pdf"
pdf(filename, width = 2.5, height = 2.5)

create.dotmap(
  x = rho,
  bg.data = pval,
  # Spot Size
  colour.scheme = c('white', 'black'),
  # Padding
  key.ylab.padding = -1.5,
  right.padding = 0.5,
  bottom.padding = 4,
  # background
  bg.alpha = 1,
  # Axes
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  yaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0,
  col.lwd = 1,
  row.lwd = 1,
  lwd = 1,
  spot.size.function = spot.size.function,
  # Legend for dots
  key = list(
    space = 'left',
    points = list(
      cex =  spot.size.function(key.sizes),
      col = default.colours(2, palette.type = 'dotmap')[2],
      pch = 19
    ),
    text = list(
      lab = as.character(key.sizes),
      cex = 0.8,
      adj = 1
    ),
    title = expression(bold(underline(rho))),
    cex.title = 0.8,
    padding.text = 6,
    background = 'white'
  ),
  colourkey = TRUE,
  at = seq(0, 10, 2), # c(10, 8, 6, 4, 2),
  colourkey.labels.at = c(0, 2, 6, 10), # seq(0, 10, 2), # c(10, 8, 6, 4, 2),
  colourkey.labels = c(
    as.expression(1),
    as.expression(changeSciNot(10^-2)),
    # as.expression(changeSciNot(10^-4)),
    as.expression(changeSciNot(10^-6)),
    # as.expression(changeSciNot(10^-8)),
    as.expression(changeSciNot(10^-10))),
  colourkey.cex = 0.8
)
dev.off()
system(paste0("cp ", filename, " ~/figures"))


# Comparison to REPIC -----------------------------------------------------

repic_gr = GRanges(seqnames = repic$V1,
                   IRanges(repic$V2, repic$V3))
repic_gr = reduce(repic_gr)

repic_overlap_helper = function(df){

  chrs = paste0("chr", c(1:22, "X", "Y"))
  df = df[df[,1] %in% chrs,]

  df_gr = GRanges(seqnames = df[,1],
                  IRanges(df[,2], df[,3]))


  nm = sum(width(intersect(repic_gr, df_gr)))
  dm = sum(width(reduce(df_gr)))
  return(c(nm/dm, nm))
}

repic_overlap = do.call(rbind, lapply(all_samples, function(this_sample){

  cat(this_sample, "\n")
  metpeak = metpeak_peaks[metpeak_peaks$sample == this_sample,]
  exomepeak = exomepeak_peaks[exomepeak_peaks$sample == this_sample,]
  metpeak_stats = repic_overlap_helper(metpeak)
  exomepeak_stats = repic_overlap_helper(exomepeak)

  data.frame(
    "Sample" = this_sample,
    "MeTPeak_p" = metpeak_stats[1],
    "ExomePeak_p" = exomepeak_stats[1],
    "MeTPeak_d" = metpeak_stats[2],
    "ExomePeak_d" = exomepeak_stats[2],
    stringsAsFactors = F
  )

}))

repic_p = melt(repic_overlap[,c("Sample", "MeTPeak_p", "ExomePeak_p")])
repic_d = melt(repic_overlap[,c("Sample", "MeTPeak_d", "ExomePeak_d")])
repic_p$variable = gsub("_.*", "", repic_p$variable)
repic_d$variable = gsub("_.*", "", repic_d$variable)
repic_df = merge(repic_p, repic_d, by = c("Sample", "variable"))

repic_df$variable = factor(repic_df$variable, levels = c("MeTPeak", "ExomePeak"))
repic_df$value.x = repic_df$value.x*100
repic_df$value.y = repic_df$value.y/10^6


filename = "figures/100_REPIC_Comparison.pdf"
pdf(filename, width = 4, height = 3.2)

create.violinplot(
  value*100 ~ variable,
  data = repic_p,
  # Axes
  ylimits = c(0, 100),
  yat = seq(0, 100, 25),
  xaxis.lab = c("exomePeak", "MeTPeak"),
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xaxis.rot = 0,
  xaxis.tck = 0,
  yaxis.tck = 0,
  ylab.label = "% Overlap with REPIC",
  ylab.cex = 1.5,
  xlab.cex = 1.5,
  xlab.label = " "
  # add.stripplot = T
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))

# create.scatterplot(
#   value.x ~ value.y,
#   data = repic_df,
#   # Point format
#   groups = repic_df$variable,
#   col = default.colours(3), # 4,
#   pch = 1,
#   alpha = 1,
#   # labels
#   ylab.label = "Overlap with peaks in the REPIC database (% bps)",
#   ylab.cex = 1,
#   xlab.label = "Coverage of Peaks (Mbps)",
#   xlab.cex = 1,
#   xaxis.cex = 1,
#   yaxis.cex = 1,
#   key = list(
#     text = list(
#       lab = c("MeTPeak (per Sample)", "exomePeak (per Sample)", "Binomial (Joint Peaks)"), # levels(repic_df$variable),
#       cex = 0.8,
#       col = 'black'
#     ),
#     points = list(
#       pch = 1,
#       col = default.colours(3), # 4
#       cex = 0.8
#     ),
#     x = 0.44,
#     y = 0.95,
#     padding.text = 2
#   )
# )
