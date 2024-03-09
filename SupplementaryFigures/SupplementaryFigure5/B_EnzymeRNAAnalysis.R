source("/cluster/home/helenzhu/code/snakemake_M6A_34_EnzymeMutations/000_HEADER.R")

library(BoutrosLab.plotting.general)
library(BoutrosLab.statistics.survival)

# Preamble ----------------------------------------------------------------
# Enzyme RNA abundance vs. Peaks and Clinical

# Loading Data ------------------------------------------------------------

# All samples
all_samples = all_samples[!all_samples %in% exclude_samples]

# Chen et al., 2019
filename = "Summary_Tables/Chen.2019.hg38/2020-10-15_RSEM_gene_TPM.txt"
chen = read.delim(filename, header = T)

# MeTPeak Peaks
filename = "Z_Matrices/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.Mask.tsv"
metpeak = read.delim(filename, header = T)
metpeak = data.frame(
  "Sample" = colnames(metpeak),
  "Peaks" = colSums(metpeak),
  stringsAsFactors = F)

# Clinical Data
clin.vars = c("sample", "pathologic_isup_grade", "summary_t", "pre_treatment_psa", "age_at_treatment", "age_at_diagnosis")
filename = "Summary_Tables/2019-08-20_500pg_SupTable1_feature_by_patient.tsv"
clinical = read.delim(filename, header = T)[,clin.vars]
clinical$age = ifelse(is.na(clinical$age_at_treatment), clinical$age_at_diagnosis, clinical$age_at_treatment)
rownames(clinical) = clinical$sample

# Grouping ISUP
clinical$pathologic_isup_grade = as.character(clinical$pathologic_isup_grade)
clinical$pathologic_isup_grade[clinical$pathologic_isup_grade %in% c("4", "5")] <- "4+"

# Grouping Clinical T Category
clinical$summary_t = substr(clinical$summary_t, 1, 2)

# IDC status
filename = "Summary_Tables/2021-05-07_m6A_idc_or_cribiform.txt"
idc = read.delim(filename, header = T)
clinical$idc = as.numeric(idc$idc_or_cribiform[match(rownames(clinical), idc$patient_id)])

# Hypoxia
filename = "Summary_Tables/Hypoxia/hypoxia_scores_buffa.txt"
hypoxia = read.delim(filename, header = T)
clinical$hypoxia = hypoxia$score[match(rownames(clinical), hypoxia$patient.IDs)]

# Survival Data
print(load("SummaryFilesNew/bcr_data.rsav"))
# [1] "bcr_data"
bcr_data$sample = gsub("[-].*", "", rownames(bcr_data))
bcr_data = unique(bcr_data[,c("bcr", "time_to_bcr", "sample")])
bcr_data = bcr_data[bcr_data$sample %in% all_samples,]
rownames(bcr_data) = bcr_data$sample
bcr_data = bcr_data[complete.cases(bcr_data),]

# Analysis ----------------------------------------------------------------

# RNA
rna = chen[enzymes$GeneID,]
rownames(rna) = enzymes$GeneName[match(rownames(rna), enzymes$GeneID)]
colnames(rna) = gsub(".F1", "", colnames(rna))
common.samples = intersect(colnames(rna), rownames(metpeak))

# MeTPeak Correlations
peak.results = do.call(rbind, lapply(1:nrow(rna), function(i){
  gene = rownames(rna)[i]
  m = cor.test(
    as.numeric(rna[gene, common.samples]),
    as.numeric(metpeak[common.samples, "Peaks"]),
    method = "spearman")
  data.frame(
    "Gene" = gene,
    "MeTPeak.p" = m$p.value,
    "MeTPeak.effectSize" = m$estimate,
    stringsAsFactors = F
  )
}))
peak.results$MeTPeak.fdr = p.adjust(peak.results$MeTPeak.p, method = "fdr")

# Clinical (Continuous, Spearman's Correlation)
clinical.continuous = function(df, clin.var, tag){
  common.samples = intersect(rownames(clinical), colnames(df))
  res = lapply(1:nrow(df), function(i){
    m = cor.test(
      as.numeric(clinical[common.samples, clin.var]),
      as.numeric(df[i, common.samples]),
      method = "spearman"
    )
    data.frame(
      "Gene" = rownames(df)[i],
      "p" = m$p.value,
      "effectSize" = m$estimate,
      stringsAsFactors = F
    )
  })
  res = do.call(rbind.data.frame, res)
  res$fdr = p.adjust(res$p, method = "fdr")
  res = res[order(res$fdr),]
  colnames(res)[colnames(res) == "p"] <- paste0(tag, ".p")
  colnames(res)[colnames(res) == "effectSize"] <- paste0(tag, ".effectSize")
  colnames(res)[colnames(res) == "fdr"] <- paste0(tag, ".fdr")
  res
}
# Age, Hypoxia, PSA
Age.results = clinical.continuous(df = rna, clin.var = "age", tag = "Age")
Hypoxia.results = clinical.continuous(df = rna, clin.var = "hypoxia", tag = "Hypoxia")
PSA.results = clinical.continuous(df = rna, clin.var = "pre_treatment_psa", tag = "PSA")

# Clinical (Categorical, One way ANOVA)
clinical.categorical = function(df, clin.var, tag){
  common.samples = intersect(rownames(clinical), colnames(df))
  res = lapply(1:nrow(df), function(i){
    tmp = data.frame(
      "clin.var" = as.character(clinical[common.samples, clin.var]),
      "rna" = as.numeric(df[i, common.samples])
    )
    res.aov <- aov(rna ~ clin.var, data = tmp)
    pval = summary(res.aov)[[1]][["Pr(>F)"]][1]
    data.frame(
      "Gene" = rownames(df)[i],
      "p" = pval,
      stringsAsFactors = F
    )
  })
  res = do.call(rbind.data.frame, res)
  res$fdr = p.adjust(res$p, method = "fdr")
  res = res[order(res$fdr),]
  colnames(res)[colnames(res) == "p"] <- paste0(tag, ".p")
  colnames(res)[colnames(res) == "fdr"] <- paste0(tag, ".fdr")
  res
}

# ISUP, T Category, IDC
ISUP.results = clinical.categorical(df = rna, clin.var = "pathologic_isup_grade", tag = "ISUP")
Summary_T.results = clinical.categorical(df = rna, clin.var = "summary_t", tag = "T")
IDC.results = clinical.categorical(df = rna, clin.var = "idc", tag = "IDC")

# BCR (Median Dichotomization)
surv = function(df, bcr, time_col = "time_to_bcr", event_col = "bcr"){

  df = t(df)
  df.totest <- as.data.frame(apply(df, 2, dichotomize.dataset))
  rownames(df.totest) <- rownames(df)
  df = df.totest

  # Identify common list of samples
  c.samples <- intersect(rownames(df), rownames(bcr));
  df <- df[c.samples,]
  bcr <- bcr[c.samples,]

  # Prepare survival object
  survobj <- Surv(bcr[,time_col], bcr[,event_col]);

  # Prepare output dataframe
  HR.table <- as.data.frame(
    matrix(
      nrow = ncol(df),
      ncol = 7,
      byrow = TRUE
    ),
    stringsAsFactors = FALSE
  );

  colnames(HR.table) <- c('HR', 'lower_95_CI_HR', 'upper_95_CI_HR', 'p_value', 'num_samples', 'logrank.p', 'loci');
  rownames(HR.table) <- colnames(df);

  # Loop for coxph
  for (i in 1:ncol(df)) {

    # Grab locus name
    this.locus <- colnames(df)[i];
    group <- as.factor(df[,i]);

    fit <- fit.coxmodel(
      groups = group,
      survobj = survobj
    );

    # Get pvalues for logrank and add it to the end
    logrank.values <- logrank.analysis(survobj, group);
    logrank.p <- logrank.values$pvalue[1];

    # Dump into output table
    HR.table[this.locus,] <- c(fit, logrank.p, this.locus);
  }

  # correct for multiple testing
  HR.table$q_value <- p.adjust(HR.table$p_value, method = 'fdr');

  # log transform the HR values
  HR.table$HR = as.numeric(HR.table$HR)
  HR.table$log2HR <- log2(HR.table$HR);

  # Return results
  return(HR.table)
}
surv.results = surv(df = rna, bcr = bcr_data, time_col = "time_to_bcr", event_col = "bcr")
surv.results = surv.results[,c("loci", "HR", "p_value", "q_value")]
colnames(surv.results) = c("Gene", "BCR.effectSize", "BCR.p", "BCR.fdr")

# Saving Results ----------------------------------------------------------

rna.results = list(
  peak.results,
  surv.results,
  ISUP.results,
  IDC.results,
  Summary_T.results,
  Age.results,
  PSA.results,
  Hypoxia.results
)
rna.results = Reduce(function(df1, df2) merge(df1, df2, by = "Gene", all = TRUE), rna.results)

save(
  rna,
  rna.results,
  file = "M_M6A_Enzyme_Mutations/RNA.m6a.enzymes.clinical.data.rsav"
)

# Write an output file for Rupert
write.table(
  rna.results,
  file = "~/m6A_enzyme_rna.tsv",
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
