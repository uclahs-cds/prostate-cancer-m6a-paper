source("/cluster/home/helenzhu/code/M6A_3_PipelineV2/000_HEADER.R")

library(reshape2)

# Preamble ----------------------------------------------------------------
# Summarizing Germline SNP Calls

setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_Germline_RNAseq/")

# Loading Data ------------------------------------------------------------

samples = read.table("/cluster/home/helenzhu/code/snakemake_M6A_3_Germline/germline_samples.txt",
                     header = T, sep = "\n", stringsAsFactors = F)

todo = expand.grid(samples$x, samples$x)

read_results = function(filename){
  tmp = readLines(filename)
  this_cols = tmp[1]
  tmp = tmp[c(1,3,4)]
  tmp = do.call(rbind, lapply(tmp, function(i) {
    a = strsplit(i, split = " ")[[1]]
    a[a != ""]
  }))
  colnames(tmp) = tmp[1,]
  tmp1 = data.frame(apply(tmp[2:3,2:ncol(tmp)], 2, as.numeric), stringsAsFactors = F)
  tmp1$threshold = tmp[2:3,1]
  tmp1
}

results_compiled = do.call(rbind, lapply(1:nrow(todo), function(i){
  cat(i, "\t")
  this_summary = paste0("./G_Comparisons/", todo[i,1], "_", todo[i,2], "/summary.txt")
  if(file.exists(this_summary)){
    cat(this_summary, "\n")
    this_results = read_results(this_summary)
    this_results$sample_1 = todo[i,1]
    this_results$sample_2 = todo[i,2]
    return(this_results)
  }
  cat("\n")
  return(NULL)
}))

samples_all = read.table("/cluster/home/helenzhu/code/snakemake_M6A_3_Germline/all_germline_samples.txt",
                             header = T, sep = "\n", stringsAsFactors = F)

todo = expand.grid(samples_all$x, samples_all$x)

results_compiled_all = do.call(rbind, lapply(1:nrow(todo), function(i){
  cat(i, "\t")
  this_summary = paste0("./H_FullComparisons/", todo[i,1], "_", todo[i,2], "/summary.txt")
  if(file.exists(this_summary)){
    cat(this_summary, "\n")
    this_results = read_results(this_summary)
    this_results$sample_1 = todo[i,1]
    this_results$sample_2 = todo[i,2]
    return(this_results)
  }
  cat("\n")
  return(NULL)
}))

save(results_compiled, results_compiled_all, file = "comparison_results_compiled.rsav")

# Create a Heatmap of the Results -----------------------------------------

print(load("comparison_results_compiled.rsav"))

sample_order = c(sort(unique(as.character(results_compiled$sample_1))),
                 sort(unique(setdiff(as.character(results_compiled_all$sample_1), as.character(results_compiled$sample_1)))))

# Covariates
sample.covariate.col = ifelse(sample_order %in% results_compiled$sample_1, default.colours(12)[1], default.colours(12)[2])
sample.covariate <- list(
  rect = list(
    col = sample.covariate.col,
    fill = sample.covariate.col,
    lwd = 1.5
  )
);

combo.cov.legend <- list(
  legend = list(
    colours = default.colours(12)[1:2],
    labels = c('European','Non-European'),
    title = 'Ancestry',
    border = 'white'
  )
)

side.legend <-legend.grob(
  legends = combo.cov.legend,
  label.cex = 1.5,
  title.cex = 1.5,
  title.just = 'left',
  title.fontface = 'bold',
  size = 2
)

# Plotting Precision
heatmap_df = results_compiled_all[results_compiled_all$threshold == "None",]
heatmap_df = acast(heatmap_df, formula = sample_1 ~ sample_2, value.var = "Precision", fun.aggregate = mean, margins = F)
heatmap_df = heatmap_df[sample_order,sample_order]

precision.hm = create.heatmap(
  # main = "Precision of RNAseq SNP Calls",
  # main.cex = 1,
  heatmap_df,
  main = "Precision",
  main.cex = 2.5,
  main.x = 0.45,
  # ylab.label = "WGS",
  # ylab.cex = 1,
  # xlab.cex = 2,
  yaxis.cex = 1,
  xaxis.cex = 1,
  clustering.method = 'none',
  colour.scheme = c("white", "red"),
  colourkey.cex = 2,
  at = seq(0.1, 0.5, 0.01),
  colourkey.labels.at = seq(0, 0.5, 0.1),
  # side covariate
  covariates = sample.covariate,
  # top covariate and covariate border specification
  covariates.top = sample.covariate,
  # covariate.legend = combo.cov.legend
  )

# Sensitivity
heatmap_df = results_compiled_all[results_compiled_all$threshold == "None",]
heatmap_df = acast(heatmap_df, formula = sample_1 ~ sample_2, value.var = "Sensitivity", fun.aggregate = mean, margins = F)
heatmap_df = heatmap_df[sample_order,sample_order]


sensitivity.hm = create.heatmap(
  # main = "Sensitivity of RNAseq SNP Calls",
  # main.cex = 1,
  heatmap_df,
  main = "Sensitivity",
  main.cex = 2.5,
  main.x = 0.45,
  # ylab.label = "WGS",
  # ylab.cex = 1,
  yaxis.cex = 1,
  xaxis.cex = 1,
  clustering.method = 'none',
  colour.scheme = c("white", "red"),
  colourkey.cex = 2,
  at = seq(0, 0.2, 0.01),
  colourkey.labels.at = seq(0, 0.2, 0.05),
  colourkey.labels = c(0, 0.05, 0.1, 0.15, 0.2),
  # side covariate
  covariates = sample.covariate,
  # top covariate and covariate border specification
  covariates.top = sample.covariate,
  # covariate.legend = combo.cov.legend
)

filename = paste0("~/figures/600_Germline_SNP_Comparison.pdf")
pdf(filename, width = 12.5, height = 8)

create.multipanelplot(
  plot.objects = list(precision.hm, sensitivity.hm),
  layout.width = 2,
  layout.height = 1,
  x.spacing = 0.7,
  ylab.label = "Whole genome sequencing",
  xlab.label = "RNA sequencing",
  ylab.cex = 2.5,
  xlab.cex = 2.5,
  # legend = list(
  #   right = list(
  #     x = 0.8,
  #     y = 1,
  #     fun = side.legend
  #   )
  # ),
  left.padding = 1,
  ylab.axis.padding = c(-5, -5),
  # xlab.axis.padding = -5,
  bottom.padding = 1,
  top.padding = 1,
  right.padding = 1
)
dev.off()
