source("/cluster/home/helenzhu/code/snakemake_M6A_26_QTLs/Rscripts/000_HEADER.R")



# Preamble ----------------------------------------------------------------
# PRS analysis for the m6A appeal

# Original: using Schumacher PRS from Katie's germline PRS paper
# Current: using Schumacher, Conti and PHS290 from Nicole's pipeline:
# Documented here: https://github.com/uclahs-cds/project-disease-ProstateTumor-PRAD-000117-GermlineBiomarkers/blob/main/script/polygenic-risk/calculate-prs-output.md

# Loading Data ------------------------------------------------------------

# Schumacher PRS - from Katie's paper
# schumacher <- read.csv("PRS/Houlahan_PRS_STable1.csv", header = T)
# schumacher <- schumacher[,c("sample", "PRS")]
# schumacher <- schumacher[schumacher$sample %in% all.euro.samples,]
# schumacher <- schumacher$PRS[match(all.euro.samples, schumacher$sample)]

# PRS data from Nicole
prs_rds <- "PRS_Nicole/2023-06-22_Houlahan_CPCG_External_UK_genotype-01_1KG-imputed_processed_prs_data.Rds"
prs_data <- readRDS(prs_rds)
prs_data <- lapply(prs_data, function(df_list) {
  df_list[['per.sample.prs.data']][,"Sample"] <- gsub("[_].*", "", rownames(df_list[['per.sample.prs.data']]))
  return(df_list)
})

# Three datasets - from each, we want the `per.sample.prs.data` and the `prs.weighted.sum` column of that data frame
prs_source <- names(prs_data)
# [1] "schumacher.prs"        "conti.multiethnic.prs" "PHS290"

# m6A Abundance Data
filename <- "Matrices/PeakCounts/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.PeakCounts.Adjusted.tsv"
m6a_adjusted_counts <- read.delim(filename, header = T)
nzeros <- apply(m6a_adjusted_counts, 1, function(x) sum(x == 0))
m6a_adjusted_counts <- m6a_adjusted_counts[nzeros == 0,]
m6a_adjusted_counts <- m6a_adjusted_counts[,all.euro.samples]

# m6A Peak Data
filename <- "Matrices/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.Mask.tsv"
m6a_binary_peaks <- read.delim(filename, header = T)
npeaks <- apply(m6a_binary_peaks, 1, function(x) sum(x))
m6a_binary_peaks <- m6a_binary_peaks[npeaks > 5,]
m6a_binary_peaks <- m6a_binary_peaks[,all.euro.samples]

# Generating a Supplementary Table ----------------------------------------

prs_list <- lapply(prs_data, function(i) i$per.sample.prs.data[,c("prs.weighted.sum", "Sample")])
prs_stable <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Sample", all.x = TRUE), prs_list)
colnames(prs_stable) <- c("Sample", "Schumacher", "Conti", "PHS290")
prs_stable <- prs_stable[prs_stable$Sample %in% all.euro.samples,]

write.table(
  prs_stable,
  file = "~/figures/prs_stable.tsv",
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

# Analysis ----------------------------------------------------------------

# Simple Spearman's correlation
simple_spearman <- function(counts, prs){

  # Filtering by variance
  var_vec <- apply(counts, 1, var)
  thresh <- summary(var_vec)["3rd Qu."]
  counts <- counts[var_vec > thresh,]

  # Spearman's Correlation
  results <- apply(counts, 1, function(x){
    res <- cor.test(x, prs, method = "pearson")
    data.frame("rho" = res$estimate, "pvalue" = res$p.value)
  })
  results <- do.call(rbind.data.frame, results)

  # Multiple testing correction
  results <- results[order(results$pvalue),]
  results$qvalue <- p.adjust(results$pvalue, method = "fdr")

  return(results)
}

# T-test
simple_ttest <- function(peaks, prs){

  # Filtering by number of samples
  n_samples <- rowSums(peaks)
  peaks <- peaks[n_samples <= 0.9*ncol(peaks) & n_samples >= 0.1*ncol(peaks),]

  # Spearman's Correlation
  results <- apply(peaks, 1, function(x){
    res <- t.test(prs[x == 1], prs[x == 0])
    data.frame(
      "n_samples" = sum(x),
      "effect_size" = mean(prs[x == 1]) - mean(prs[x == 0]),
      "pvalue" = res$p.value
    )
  })
  results <- do.call(rbind.data.frame, results)

  # Multiple testing correction
  results <- results[order(results$pvalue),]
  results$qvalue <- p.adjust(results$pvalue, method = "fdr")

  return(results)
}

peak_results <- lapply(prs_data, function(df_list){
  this_prs <- df_list$per.sample.prs.data$prs.weighted.sum[match(all.euro.samples, df_list$per.sample.prs.data$Sample)]
  simple_ttest(peaks = m6a_binary_peaks, prs = this_prs)
})

abundance_results <- lapply(prs_data, function(df_list){
  this_prs <- df_list$per.sample.prs.data$prs.weighted.sum[match(all.euro.samples, df_list$per.sample.prs.data$Sample)]
  simple_spearman(counts = m6a_adjusted_counts, prs = this_prs)
})

# Significant results stable ----------------------------------------------

sig_conti_results <- peak_results[["conti.multiethnic.prs"]]
sig_conti_results <- sig_conti_results[sig_conti_results$qvalue < 0.1,]
write.table(
  sig_conti_results,
  file = "~/figures/conti_sig_stable.tsv",
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)

# Plotting ----------------------------------------------------------------

# Generating labels
generate_labels <- function(ids){
  split_ids <- strsplit(ids, split = "[:]")
  sapply(split_ids, function(spl_id) paste0(gene_conversion[spl_id[1]], " | ", spl_id[2]))
}

# Histogram of P-values
pval_peak_histogram <- function(peak_res, prs_name) {
  plt <- create.histogram(
    peak_res$pvalue,
    # type = "count",
    breaks = seq(0, 1, 0.1),
    # Axes
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xlab.cex = 1.5,
    ylab.cex = 1.5,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xlab.label = paste0(prs_name, " P value"),
    ylab.label = "Percent",
    xlimits = c(0, 1),
    ylimits = c(0, 15),
    xat = seq(0, 1, 0.2),
    yat = seq(0, 15, 5)
    # main = prs_name,
    # main.cex = 1.5
  )
  return(plt)
}

pval_abundance_histogram <- function(abundance_res, prs_name){
  plt <- create.histogram(
    abundance_res$pvalue,
    # type = "count",
    breaks = seq(0, 1, 0.1),
    # Axes
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xlab.cex = 1.5,
    ylab.cex = 1.5,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xlab.label = paste0(prs_name, " P value"),
    ylab.label = "Percent",
    xlimits = c(0, 1),
    ylimits = c(0, 20),
    xat = seq(0, 1, 0.2),
    yat = seq(0, 20, 5)
    # main = prs_name,
    # main.cex = 1.5
  )
  return(plt)
}

# Volcano plots
volcano_peak <- function(peak_res, prs_name, label_sig_points = T){

  sig_points <- peak_res[peak_res$qvalue < 0.1,]
  sig_points$labels <- generate_labels(rownames(sig_points))

  plt <- create.scatterplot(
    -log10(pvalue) ~ effect_size,
    peak_res,
    # Axes
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xlab.cex = 1.5,
    ylab.cex = 1.5,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xlab.label = bquote(bold(Delta ~ .(prs_name) ~"PRS")),
    ylab.label = expression(bold("-log"[10]*" P value")),
    xlimits = c(-1, 1),
    ylimits = c(0, 6.5),
    xat = seq(-1, 1, 0.5),
    yat = seq(0, 6, 2),
    # main = prs_name,
    # main.cex = 1.5,
    # Adding sig points
    add.points = TRUE,
    points.x = sig_points$effect_size,
    points.y = -log10(sig_points$pvalue),
    points.pch = 19,
    points.col = 'blue',
    points.col.border = 'blue',
    points.cex = 1,
    # Adding labels
    add.text = TRUE,
    text.labels = sig_points$labels,
    text.x = sig_points$effect_size,
    text.y = -log10(sig_points$pvalue) + 0.25,
    text.col = 'blue',
    text.cex = 1.2,
    text.fontface = 'bold'
    # text.guess.labels = FALSE,
    # text.guess.skip.labels = TRUE,
    # text.guess.ignore.radius = FALSE,
    # text.guess.ignore.rectangle = FALSE,
    # text.guess.radius.factor = 1,
    # text.guess.buffer.factor = 1,
    # text.guess.label.position = 180
  )

  return(plt)
}

volcano_abundance <- function(abundance_res, prs_name, label_sig_points = T){

  sig_points <- abundance_res[abundance_res$qvalue < 0.1,]
  sig_points$labels <- generate_labels(rownames(sig_points))

  plt <- create.scatterplot(
    -log10(pvalue) ~ rho,
    abundance_res,
    # Axes
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xlab.cex = 1.5,
    ylab.cex = 1.5,
    xaxis.tck = 0,
    yaxis.tck = 0,
    xlab.label = bquote(bold(.(prs_name) ~"PRS" ~ rho)),
    ylab.label = expression(bold("-log"[10]*" P value")),
    xlimits = c(-0.5, 0.5),
    ylimits = c(0, 6.5),
    xat = seq(-0.5, 0.5, 0.25),
    yat = seq(0, 6, 2),
    # main = prs_name,
    # main.cex = 1.5,
    # Adding sig points
    add.points = TRUE,
    points.x = sig_points$rho,
    points.y = -log10(sig_points$pvalue),
    points.pch = 19,
    points.col = 'blue',
    points.col.border = 'blue',
    points.cex = 1,
    # Adding labels
    add.text = TRUE,
    text.labels = sig_points$labels,
    text.x = sig_points$rho,
    text.y = -log10(sig_points$pvalue) + 0.25,
    text.col = 'blue',
    text.cex = 1.2,
    text.fontface = 'bold',
    # text.guess.labels = FALSE,
    # text.guess.skip.labels = TRUE,
    # text.guess.ignore.radius = FALSE,
    # text.guess.ignore.rectangle = FALSE,
    # text.guess.radius.factor = 1,
    # text.guess.buffer.factor = 1,
    # text.guess.label.position = NULL,
  )
  return(plt)
}

# Generating Plots --------------------------------------------------------

filename = "~/figures/500_PRS_m6A_peak_ttest.pdf"
pdf(filename, width = 10, height = 8)

plot_objects <- list(
  pval_peak_histogram(peak_results$schumacher.prs, prs_name = "Schumacher"),
  pval_peak_histogram(peak_results$conti.multiethnic.prs, prs_name = "Conti"),
  pval_peak_histogram(peak_results$PHS290, prs_name = "PHS290"),
  volcano_peak(peak_results$schumacher.prs, prs_name = "Schumacher"),
  volcano_peak(peak_results$conti.multiethnic.prs, prs_name = "Conti"),
  volcano_peak(peak_results$PHS290, prs_name = "PHS290")
)

create.multipanelplot(
  plot_objects,
  # Layout
  layout.height = 2,
  layout.width = 3
)

dev.off()

filename = "~/figures/500_PRS_m6A_abundance_spearman.pdf"
pdf(filename, width = 14, height = 10)

plot_objects <- list(
  pval_abundance_histogram(abundance_results$schumacher.prs, prs_name = "Schumacher"),
  pval_abundance_histogram(abundance_results$conti.multiethnic.prs, prs_name = "Conti"),
  pval_abundance_histogram(abundance_results$PHS290, prs_name = "PHS290"),
  volcano_abundance(abundance_results$schumacher.prs, prs_name = "Schumacher"),
  volcano_abundance(abundance_results$conti.multiethnic.prs, prs_name = "Conti"),
  volcano_abundance(abundance_results$PHS290, prs_name = "PHS290")
)

create.multipanelplot(
  plot_objects,
  # Layout
  layout.height = 2,
  layout.width = 3
)

dev.off()

filename = "~/figures/500_PRS_m6A_peak_ttest.individual.pdf"
pdf(filename, width = 5, height = 5)

pval_peak_histogram(peak_results$schumacher.prs, prs_name = "Schumacher")
pval_peak_histogram(peak_results$conti.multiethnic.prs, prs_name = "Conti")
pval_peak_histogram(peak_results$PHS290, prs_name = "PHS290")
volcano_peak(peak_results$schumacher.prs, prs_name = "Schumacher")
volcano_peak(peak_results$conti.multiethnic.prs, prs_name = "Conti")
volcano_peak(peak_results$PHS290, prs_name = "PHS290")

dev.off()

# Paper plots -------------------------------------------------------------

sig_results <- peak_results[["conti.multiethnic.prs"]]
sig_results <- sig_results[sig_results$qvalue < 0.1,]
plotting_data <- data.frame(
  "Peak" = as.numeric(m6a_binary_peaks[rownames(sig_results),all.euro.samples]),
  "Conti" = prs_data[["conti.multiethnic.prs"]]$per.sample.prs.data$prs.weighted.sum[match(all.euro.samples, prs_data[["conti.multiethnic.prs"]]$per.sample.prs.data$Sample)]
)
plotting_data$Peak <- ifelse(plotting_data$Peak == 1, "Peak", "No Peak")
plotting_data$Peak <- factor(plotting_data$Peak, levels = c("No Peak", "Peak"))

this_lab <- rownames(sig_results)[1]
this_lab <- strsplit(this_lab, split = "[::]")[[1]]
this_lab <- paste0(gene_conversion[this_lab[[1]]], " | ", this_lab[[2]])

qval <- formatC(sig_results$qvalue[1], digits = 2)
delta <- formatC(sig_results$effect_size[1], digits = 2)

filename ="~/figures/500_PRS_sig_results_boxplot.pdf"
pdf(filename, width = 5, height = 6)

create.boxplot(
  Conti ~ Peak,
  plotting_data,
  # Stripplot
  add.stripplot = T,
  # Formatting
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xlab.cex = 1.5,
  xlab.label = this_lab,
  ylab.cex = 1.5,
  ylab.label = "Conti PRS",
  # Title
  # main.cex = 1,
  # main = this_lab,
  # Adding P value
  add.text = TRUE,
  text.labels = as.expression(
    bquote(
      atop(
      "Q < " ~ .(qval),
      Delta ~ "PRS =" ~ .(delta)
  ))),
  text.x = 2,
  text.y = 24.5,
  text.anchor = 'centre',
  text.col = 'black',
  text.cex = 1.7,
  text.fontface = 1,
)

dev.off()
