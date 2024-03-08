source("/cluster/home/helenzhu/code/snakemake_M6A_11_SubsamplingBams/RScripts/000_HEADER.R")

# Preamble ----------------------------------------------------------------
# Collecting data for different Samples & Seeds for MeTPeak, exomePeak
# 2020-09-15: Excluding MACS2 Results

# Loading Data ------------------------------------------------------------

# Defining Parameteres
samples = c("CPCG0464", "CPCG0562", "CPCG0382")
seeds = c('314159', '271828', '1618')
percent = c("0.5", "0.25", "0.1", "0.01")
peak_methods = c("exomePeak", "MeTPeak") # MACS2

# Defining Files
todo = expand.grid(seeds, percent, samples)
files = paste0(todo[,3], "_IP_", todo[,2], "_", todo[,1])

# Exomepeak & MeTPeak
exomepeak_metpeak_files_bed6 = paste0(files, "/peak.bed6")
exomepeak_metpeak_oracle_bed6 = paste0(samples, "/peak.bed6")
exomepeak_metpeak_files_xls = paste0(files, "/peak.xls")
exomepeak_metpeak_oracle_xls = paste0(samples, "/peak.xls")
exomepeak_dir = "O_PeaksSim/exomePeak_Peaks/"
metpeak_dir = "O_PeaksSim/MeTPeak_Peaks/"
exomepeak_oracle_dir = "G_exomePeak_Peaks_SS/"
metpeak_oracle_dir = "F_MeTPeak_Peaks_SS/"

# Loading Metpeak Data
load_peaks = function(filename, directory, peak_method, oracle = F, xls = F){

  # Print filename
  cat(filename, "\n")

  # Importing file
  tmp = read.table(paste0(directory, filename),
                   sep = "\t",
                   header = xls,
                   stringsAsFactors = F,
                   comment.char = "#")

  # Adding columns identifiers
  tag = strsplit(filename, split = "[_]|[/]")[[1]]

  if(oracle){
    tmp$Sample = tag[1]
    tmp$Experiment = "IP"
    tmp$Proportion = "1"
    tmp$Seed = "0"
  } else{
    tmp$Sample = tag[1]
    tmp$Experiment = tag[2]
    tmp$Proportion = tag[3]
    tmp$Seed = tag[4]
  }

  tmp$Method = peak_method
  tmp
}

# Loading ExomePeak
exomepeak_b6 = do.call(rbind, lapply(exomepeak_metpeak_files_bed6, function(i)
  load_peaks(filename = i, directory = exomepeak_dir, peak_method = "exomePeak")))
exomepeak_oracle_b6 = do.call(rbind, lapply(exomepeak_metpeak_oracle_bed6, function(i)
  load_peaks(filename = i, directory = exomepeak_oracle_dir, peak_method = "exomePeak", oracle = T)))

exomepeak_xls = do.call(rbind, lapply(exomepeak_metpeak_files_xls, function(i)
  load_peaks(filename = i, directory = exomepeak_dir, peak_method = "exomePeak", xls = T)))
exomepeak_oracle_xls = do.call(rbind, lapply(exomepeak_metpeak_oracle_xls, function(i)
  load_peaks(filename = i, directory = exomepeak_oracle_dir, peak_method = "exomePeak", oracle = T, xls = T)))

# Loading MeTPeak
metpeak_b6 = do.call(rbind, lapply(exomepeak_metpeak_files_bed6, function(i)
  load_peaks(filename = i, directory = metpeak_dir, peak_method = "MeTPeak")))
metpeak_oracle_b6 = do.call(rbind, lapply(exomepeak_metpeak_oracle_bed6, function(i)
  load_peaks(filename = i, directory = metpeak_oracle_dir, peak_method = "MeTPeak", oracle = T)))

metpeak_xls = do.call(rbind, lapply(exomepeak_metpeak_files_xls, function(i)
  load_peaks(filename = i, directory = metpeak_dir, peak_method = "MeTPeak", xls = T)))
metpeak_oracle_xls = do.call(rbind, lapply(exomepeak_metpeak_oracle_xls, function(i)
  load_peaks(filename = i, directory = metpeak_oracle_dir, peak_method = "MeTPeak", oracle = T, xls = T)))

# Compiling Peaks
peaks_compiled_b6 = rbind(exomepeak_b6, exomepeak_oracle_b6, metpeak_b6, metpeak_oracle_b6)
peaks_compiled_xls = rbind(exomepeak_xls, exomepeak_oracle_xls, metpeak_xls, metpeak_oracle_xls)

# Number of reads
print(load("Summary_Tables/hg38.STAR.output.rsav"))
# [1] "star_output"
star_output = star_output[star_output$Sample %in% samples, c("Sample", "UniquelyMapped.Count.IP", "MultiMapped.Count.IP")]
star_output$total = star_output$UniquelyMapped.Count.IP + star_output$MultiMapped.Count.IP
star_output = star_output[,c("Sample", "total")]

# Making Some Nice Plots --------------------------------------------------

# 1. Number of Peaks vs. Percent/Number of Reads
create_peaks.vs.reads = function(df, this_method, this_main, this_ylab, this_xlab, read_count = F){

  df = df[df$Method == this_method,]
  df$Proportion = as.numeric(df$Proportion)
  df$count = 1
  df = aggregate(count ~ Sample + Proportion + Seed + Method, df, FUN = sum)
  df$Sample = factor(df$Sample, levels = samples)

  if(read_count){
    df = merge(df, star_output, by = "Sample")
    df$Proportion = df$Proportion*df$total/10^6
  }

  #  this_main.x = ifelse(this_main == "MeTPeak", 0.26, 0.3)

  plt = create.scatterplot(
    count/1000 ~ Proportion,
    data = df,
    groups = df$Sample,
    # Axis labels
    xaxis.cex = 0.8,
    yaxis.cex = 0.8,
    xlab.cex = 0.8,
    ylab.cex = 0.8,
    ylab.label = this_ylab,
    xlab.label = this_xlab,
    # ylimits = c(0, 15),
    # yat = seq(0, 15, 5),
    # Main
    # main = this_main,
    # main.cex = 0.5,
    # main.x = this_main.x,
    # main.y = -5,
    # main.just = 'left',
    # Tcks
    xaxis.tck = 0,
    yaxis.tck = 0,
    # Lines and pch's
    pch = 1,
    type =  c("a", "p"),
    # Colour & Legend
    col = default.colours(3)
  )
  plt
}

# 2. Percent overlap with oracle & correlation of results
overlap_oracle = function(id, proportion, seed, method){

  # id = "CPCG0464"
  # proportion = 0.5
  # seed = 314159
  # method = "exomePeak"

  cat(id, "\t", proportion, "\t", seed, "\t", method, "\n")
  peaks_oracle = peaks_compiled_b6[peaks_compiled_b6$Method == method &
                                   peaks_compiled_b6$Sample == id,]
  peaks_sample = peaks_compiled_b6[peaks_compiled_b6$Method == method &
                                   peaks_compiled_b6$Sample == id &
                                   peaks_compiled_b6$Seed == seed &
                                   peaks_compiled_b6$Proportion == proportion,]

  oracle_gr = GRanges(seqnames = peaks_oracle$V1,
                      IRanges(peaks_oracle$V2, peaks_oracle$V3))
  peaks_gr = GRanges(seqnames = peaks_sample$V1,
                     IRanges(peaks_sample$V2, peaks_sample$V3))

  # Overlap
  this_overlap = sum(width(reduce(intersect(oracle_gr, peaks_gr))))
  this_oracle = sum(width(reduce(oracle_gr)))
  this_peaks = sum(width(reduce(peaks_gr)))
  percent_oracle = this_overlap/this_oracle
  percent_peaks = this_overlap/this_peaks

  # Results
  data.frame(
    "Sample" = id,
    "Proportion" = proportion,
    "Seed" = seed,
    "Method" = method,
    "Percent_Oracle" = percent_oracle,
    "Percent_Peaks" = percent_peaks,
    stringsAsFactors = F
  )
}

correlation_oracle = function(id, proportion, seed, method){

  # id = "CPCG0464"
  # proportion = 0.5
  # seed = 314159
  # method = "exomePeak"

  cat(id, "\t", proportion, "\t", seed, "\t", method, "\n")
  peaks_oracle = peaks_compiled_xls[peaks_compiled_xls$Method == method &
                                     peaks_compiled_xls$Sample == id,]
  peaks_sample = peaks_compiled_xls[peaks_compiled_xls$Method == method &
                                     peaks_compiled_xls$Sample == id &
                                     peaks_compiled_xls$Seed == seed &
                                     peaks_compiled_xls$Proportion == proportion,]

  oracle_gr = GRanges(seqnames = peaks_oracle$chr,
                      IRanges(peaks_oracle$chromStart, peaks_oracle$chromEnd))
  peaks_gr = GRanges(seqnames = peaks_sample$chr,
                     IRanges(peaks_sample$chromStart, peaks_sample$chromEnd))

  # Correlation
  p_oracle = peaks_oracle[queryHits(findOverlaps(oracle_gr, peaks_gr)), c("score", "fold_enrchment")]
  p_sample = peaks_sample[subjectHits(findOverlaps(oracle_gr, peaks_gr)), c("score", "fold_enrchment")]
  p.spearman = cor(p_oracle$score, p_sample$score, method = "spearman")
  fc.spearman = cor(p_oracle$fold_enrchment, p_sample$fold_enrchment, method = "spearman")

  # Results
  data.frame(
    "Sample" = id,
    "Proportion" = proportion,
    "Seed" = seed,
    "Method" = method,
    "P_Spearman" = p.spearman,
    "FC_Spearman" = fc.spearman,
    stringsAsFactors = F
  )
}

plot_stats = function(df, this_method, this_main, this_ylab, this_xlab, spearman = F, read_count = F){

  # df = overlap_peaks
  # this_method = "exomePeak"
  # this_main = ""
  # this_ylab = ""
  # this_xlab = ""
  # spearman = F
  # read_count = F

  df = df[df$Method == this_method,]

  if(spearman){
    df_oracle = data.frame(
      "Sample" = samples,
      "Proportion" = "1",
      "Seed" = "0",
      "Method" = this_method,
      "P_Spearman" = 1,
      "FC_Spearman" = 1,
      stringsAsFactors = F
    )

    df = rbind(df, df_oracle)
    df$Proportion = as.numeric(df$Proportion)

    df$Y = df$FC_Spearman
    this_yat = seq(0, 1, 0.2)
    this_ylimits = c(0, 1.05)

  } else{
    df_oracle = data.frame(
      "Sample" = samples,
      "Proportion" = "1",
      "Seed" = "0",
      "Method" = this_method,
      "Percent_Oracle" = 1,
      "Percent_Peaks" = 1,
      stringsAsFactors = F
    )

    df = rbind(df, df_oracle)
    df$Proportion = as.numeric(df$Proportion)

    df$Y = df$Percent_Oracle*100
    this_yat = seq(0, 100, 20)
    this_ylimits = c(0, 105)
  }

  if(read_count){
    df = merge(df, star_output, by = "Sample")
    df$Proportion = df$Proportion*df$total/10^6
  }

  this_main.x = ifelse(this_main == "MeTPeak", 0.26, 0.3)

  plt = create.scatterplot(
    Y ~ Proportion,
    data = df,
    groups = df$Sample,
    # Axis labels
    xaxis.cex = 0.5,
    yaxis.cex = 0.5,
    xlab.cex = 0.5,
    ylab.cex = 0.5,
    ylab.label = this_ylab,
    xlab.label = this_xlab,
    ylimits = this_ylimits,
    yat = this_yat,
    # Main
    main = this_main,
    main.cex = 0.5,
    main.x = this_main.x,
    main.y = -5,
    main.just = 'left',
    # Tcks
    xaxis.tck = 0,
    yaxis.tck = 0,
    # Lines and pch's
    pch = 1,
    type =  c("a", "p"),
    # Colour & Legend
    col = default.colours(3)
  )
  plt
}

todo = expand.grid(seeds, percent, samples, peak_methods, stringsAsFactors = F)
overlap_peaks = do.call(
  rbind, lapply(1:nrow(todo),
                function(i) overlap_oracle(id = todo[i, "Var3"],
                                           proportion = todo[i, "Var2"],
                                           seed = todo[i, "Var1"],
                                           method = todo[i, "Var4"])))
correlation_peaks = do.call(
  rbind, lapply(1:nrow(todo),
                function(i) correlation_oracle(id = todo[i, "Var3"],
                                               proportion = todo[i, "Var2"],
                                               seed = todo[i, "Var1"],
                                               method = todo[i, "Var4"])))
# Compiling into a full figure --------------------------------------------

return_plot = function(this_method, plottype){

  cat(this_method, "\t", plottype, "\n")
  if(plottype == 1){
    this_ylab = ifelse(this_method == "exomePeak",  expression(bold("# exomePeak peaks (per 10"^"3"*")")), expression(bold("# MeTPeak peaks (per 10"^"3"*")")))
    create_peaks.vs.reads(df = peaks_compiled_xls,
                          this_method = this_method,
                          this_main = this_method,
                          this_ylab = this_ylab,
                          this_xlab = "Proportion of Total Reads")
  } else if (plottype == 2){
    this_ylab = ifelse(this_method == "exomePeak", expression(bold("exomePeak peaks (10"^"3"*")")), expression(bold("MeTPeak peaks (10"^"3"*")")))
    this_xlab = expression(bold("Reads (10"^"6"*")"))
    create_peaks.vs.reads(df = peaks_compiled_xls,
                          this_method = this_method,
                          this_main = this_method,
                          this_ylab = this_ylab,
                          this_xlab = this_xlab,
                          read_count = T)
  } else if (plottype == 3){
    this_ylab = ifelse(this_method == "exomePeak", "% Overlap", "")
    plot_stats(df = overlap_peaks,
               this_method = this_method,
               this_main = "",
               this_ylab = this_ylab,
               this_xlab = "Proportion of Reads")
  } else if (plottype == 4){
    this_ylab = ifelse(this_method == "exomePeak", "% Overlap", "")
    plot_stats(df = overlap_peaks,
               this_method = this_method,
               this_main = "",
               this_ylab = this_ylab,
               this_xlab = "",
               read_count = T)
  } else if (plottype == 5){
    this_ylab = ifelse(this_method == "exomePeak", "Spearman's Rho", "")
    plot_stats(df = correlation_peaks,
               this_method = this_method,
               this_main = "",
               this_ylab = this_ylab,
               this_xlab = "",
               spearman = T)
  } else if (plottype == 6){
    this_ylab = ifelse(this_method == "exomePeak", "Spearman's Rho", "")
    plot_stats(df = correlation_peaks,
               this_method = this_method,
               this_main = "",
               this_ylab = this_ylab,
               this_xlab = "",
               spearman = T,
               read_count = T)
  }
}


# Paper Figure ------------------------------------------------------------

filename = "figures/100_PeaksSim.Paper.pdf"
pdf(filename, width = 6, height = 2.5)

# Legend
this_legend <- list(
  legend = list(
    colours = default.colours(3),
    labels = samples,
    border = 'white',
    title = 'Samples',
    pch = 1
  )
)

legend.grob <- legend.grob(
  legends = this_legend,
  title.just = 'left',
  label.cex = 0.7,
  title.cex = 0.7
)

plot_types = c(2)
plot_todo = expand.grid(peak_methods, plot_types, stringsAsFactors = F)
plotting.objects = lapply(1:nrow(plot_todo), function(i)
  return_plot(this_method = plot_todo[i, "Var1"], plottype = plot_todo[i, "Var2"]))

create.multipanelplot(
  plot.objects = plotting.objects,
  layout.width = 2,
  layout.height = 1,
  x.spacing = 0,
  y.spacing = -2,
  xlab.cex = 0,
  ylab.cex = 0,
  # main.x = '',
  # main.y = '',
  main.cex = 0,
  plot.objects.heights = 1,
  plot.objects.widths = c(1, 1),
  legend = list(
    right = list(
      x = 0.10,
      y = 0.50,
      fun = legend.grob
    )
  ),
  left.legend.padding = 0,
  right.legend.padding = 0,
  bottom.legend.padding = 0,
  top.legend.padding = 0
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))

# plot_types = c(1)
# plot_todo = expand.grid(peak_methods, plot_types, stringsAsFactors = F)
# plotting.objects = lapply(1:nrow(plot_todo), function(i)
#   return_plot(this_method = plot_todo[i, "Var1"], plottype = plot_todo[i, "Var2"]))
#
# create.multipanelplot(
#   plot.objects = plotting.objects,
#   layout.width = 1,
#   layout.height = 2,
#   x.spacing = -0.5,
#   y.spacing = -2,
#   xlab.cex = 0,
#   ylab.cex = 0,
#   # main.x = '',
#   # main.y = '',
#   main.cex = 0,
#   plot.objects.heights = c(1.0, 1.1),
#   plot.objects.widths  = 1,
#   legend = list(
#     right = list(
#       x = 0.10,
#       y = 0.50,
#       fun = legend.grob
#     )
#   )
# )


# Compiling Figure --------------------------------------------------------

filename = "figures/100_PeaksSim.pdf"
pdf(filename, width = 5, height = 5)

# Legend
this_legend <- list(
  legend = list(
    colours = default.colours(3),
    labels = samples,
    border = 'white',
    title = 'Samples',
    pch = 1
  )
)

legend.grob <- legend.grob(
  legends = this_legend,
  title.just = 'left',
  label.cex = 0.5,
  title.cex = 0.5
)

plot_types = seq(2, 6, 2)
plot_todo = expand.grid(peak_methods, plot_types, stringsAsFactors = F)
plotting.objects = lapply(1:nrow(plot_todo), function(i)
  return_plot(this_method = plot_todo[i, "Var1"], plottype = plot_todo[i, "Var2"]))

create.multipanelplot(
  plot.objects = plotting.objects,
  layout.width = 2,
  layout.height = 3,
  x.spacing = -0.5,
  y.spacing = -2,
  xlab.label = "Number of Reads (per Million)",
  xlab.cex = 0.8,
  ylab.label = "Metric",
  ylab.cex = 0.8,
  main = 'Testing the Effects of Downsampling BAMs on Peak Calling',
  main.cex = 0.8,
  # main.x = '',
  # main.y = '',
  plot.objects.heights = c(1.05, 1, 1),
  plot.objects.widths = c(1.07, 1),
  legend = list(
    right = list(
      x = 0.10,
      y = 0.50,
      fun = legend.grob
    )
  )
)

plot_types = seq(1, 5, 2)
plot_todo = expand.grid(peak_methods, plot_types, stringsAsFactors = F)
plotting.objects = lapply(1:nrow(plot_todo), function(i)
  return_plot(this_method = plot_todo[i, "Var1"], plottype = plot_todo[i, "Var2"]))

create.multipanelplot(
  plot.objects = plotting.objects,
  layout.width = 2,
  layout.height = 3,
  x.spacing = -0.5,
  y.spacing = -2,
  xlab.label = "Proportion of Total Reads",
  xlab.cex = 0.8,
  ylab.label = "Metric",
  ylab.cex = 0.8,
  main = 'Testing the Effects of Downsampling BAMs on Peak Calling',
  main.cex = 0.8,
  # main.x = '',
  # main.y = '',
  plot.objects.heights = c(1.05, 1, 1),
  plot.objects.widths = c(1.07, 1),
  legend = list(
    right = list(
      x = 0.10,
      y = 0.50,
      fun = legend.grob
    )
  )
)

dev.off()
system(paste0("cp ", filename, " ~/figures"))
