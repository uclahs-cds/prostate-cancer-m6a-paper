source("000_HEADER.R")

library(peer)

# Preamble ----------------------------------------------------------------
# Trying out the software PEER


# Samples -----------------------------------------------------------------

# All European Samples
all.euro.samples = readLines("all.euro.samples.tsv")

# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly = T)
meth = args[1]
tag = args[2]
n_components = strtoi(args[3])

# meth = "MeTPeak"
# tag = 'V9_MetricVoting_OptimizationFunction_MaxGapOnly'
# n_components = 1

# Loading Expression Data
expr.file = paste0("Input/", meth, ".", tag, ".Expression.tsv")
expr = read.table(
  file = expr.file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

# Formatting
rownames(expr) = expr$peak
expr = expr[,all.euro.samples]
expr = as.matrix(expr)
expr = apply(expr, 2, as.numeric)
expr = t(expr)

# Loading Covariates - Age
age.file = "Summary_Tables/2019-08-20_500pg_SupTable1_feature_by_patient.tsv"
age = read.delim(age.file, header = T)
age$age = ifelse(is.na(age$age_at_treatment), age$age_at_diagnosis, age$age_at_treatment)
rownames(age) = age$patient_id
age = age[all.euro.samples, "age", drop = F]
age = scale(age, center = T, scale = T)

# Loading Covariates - Plink Genotype
plink.file = "Plink/CPCGENE_filtered.covariates"
plink = read.delim(plink.file, header = T)
plink = t(plink)

covs = cbind(age, plink)

# Zero Peer Covariates ----------------------------------------------------

if(n_components == 0){
  filename = paste0("covariates/", meth , ".", tag, "_", n_components, "_peer.tsv")
  colnames(covs) = paste0("Cov", 1:ncol(covs))
  write.table(
    t(covs),
    file = filename,
    sep = "\t",
    col.names = T,
    row.names = T,
    quote = F
  )
} else {

  # Setting Up PEER ---------------------------------------------------------

  # Initializing Model
  model = PEER()

  # Setting the number of components
  PEER_setNk(model, n_components)

  # Checking for the number of hidden components
  PEER_getNk(model)

  # Adding Expression Data
  PEER_setPhenoMean(model, expr)

  # Checking for the dimension of the matrix (Sample, row X Gene, column)
  dim(PEER_getPhenoMean(model))

  # Adding Expression Means as a Covariate
  # PEER_setAdd_mean(model, TRUE)

  # Setting Covariates
  PEER_setCovariates(model, as.matrix(covs))

  # Additional Parameters ---------------------------------------------------

  # Number of iterations
  # PEER_setNmax_iterations(model, 100)

  # Threshold for increase in variational lower bound
  # PEER_setTolerance(model, 1)

  # Threshold for change in variance
  # PEER_setVarTolerance(model, 0.1)

  # Initializing priors for gamma for noise and precision, a parameter, b parameter
  # PEER_setPriorAlpha(model,0.001,0.1)
  # PEER_setPriorEps(model,0.1,10.)

  # Running the Model -------------------------------------------------------

  # Running the model
  PEER_update(model)

  # Extracting Results
  factors = PEER_getX(model)
  dim(factors)

  weights = PEER_getW(model)
  dim(weights)

  precision = PEER_getAlpha(model)
  dim(precision)

  residuals = PEER_getResiduals(model)
  dim(residuals)

  # plot(precision)

  # Saving the Results ------------------------------------------------------

  # factors
  filename = paste0("covariates/", meth , ".", tag, "_", n_components, "_peer.tsv")
  colnames(factors) = paste0("Cov", 1:ncol(factors))
  rownames(factors) = all.euro.samples
  write.table(
    t(factors),
    file = filename,
    sep = "\t",
    col.names = T,
    row.names = T,
    quote = F
  )

  # weights
  rownames(weights) = colnames(expr)
  colnames(weights) = colnames(factors)

  # precision
  rownames(precision) = colnames(factors)
  colnames(precision) = "Precision"

  # residuals
  rownames(residuals) = all.euro.samples
  colnames(residuals) = colnames(expr)

  # Saving file
  filename = paste0("covariates/", meth , ".", tag, "_", n_components, "_peer.rsav")
  save(factors,
       weights,
       precision,
       residuals,
       file = filename)

}
