source("000_HEADER.R")
source("RADAR_diffIP_parallel_handle_NA.R")

library(RADAR)

# Preamble ----------------------------------------------------------------
# This does a differential methylation analysis on
# 1. Tumour-Normal Peaks
# 2. Xenografts

# Stealing code from Rupert for this one:
# https://github.com/uclahs-cds/project-ProstateCancer-m6A/blob/rhughwhite/somatic_driver_events/01_driver_m6A_enzyme_differential_testing.R

# Loading Data ------------------------------------------------------------

# Tumour Normal
tumour.normal = read.delim("H_PeakCounts/TumourNormal.Adjusted.Counts.tsv", header = T)

# Xenografts
xenografts = read.delim( "H_PeakCounts/Xenograft.Adjusted.Counts.tsv", header = T)

# RADAR Test --------------------------------------------------------------

test.mutation.RADAR <- function(
  mutation,
  mutation.load = NULL,
  radar.obj,
  threads
) {
  if (!is.null(mutation.load)) {
    variable(radar.obj) <- data.frame(mutation, mutation.load);
  } else {
    variable(radar.obj) <- data.frame(mutation);
  }
  radar.test <- diffIP_parallel_NA(radar.obj, thread = threads)@test.est;
  errors <- grep('Error', radar.test[, 'p_value']);
  if (length(errors) > 0 ) {
    radar.test[errors, ] <- NA;
    radar.test <- apply(radar.test, 2, as.numeric);
    row.names(radar.test) <- row.names(radar.obj@ip_adjExpr_filtered);
  }
  return(radar.test);
}

run.test = function(
  m6A,
  covariate
){
  radar.obj <- MeRIP.RADAR();
  samplenames(radar.obj) <- colnames(m6A);
  radar.obj@fdr.method <- 'qvalue';
  radar.obj@ip_adjExpr_filtered <- as.matrix(m6A);


  driver.results = test.mutation.RADAR(
    mutation = covariate,
    radar.obj = radar.obj,
    threads = 1
  )

  return(driver.results)
}

# Running the test on xenograft and tumour-normal -------------------------

tn.covariate = as.numeric(grepl("^N", colnames(tumour.normal)))
tn.results = run.test(
  m6A = tumour.normal,
  covariate = tn.covariate
)
tn.results = data.frame(tn.results)

x.covariate = c(0, 0, 1, 1)
x.results = run.test(
  m6A = xenografts,
  covariate = x.covariate
)
x.results = data.frame(x.results)

# Writing Result Tables ---------------------------------------------------

write.table(
  tn.results,
  file = "DifferentialMethylation/TumourNormal.results.tsv",
  sep = '\t',
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)

write.table(
  x.results,
  file = "DifferentialMethylation/Xenograft.results.tsv",
  sep = '\t',
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)
