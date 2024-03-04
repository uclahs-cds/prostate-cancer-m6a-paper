source("/cluster/home/helenzhu/code/snakemake_M6A_8_PipelineRebuild/Rscripts/000_HEADER.R")

library(survival)
library(BoutrosLab.plotting.general)

# Preamble ----------------------------------------------------------------

# TODO:
# Decide if we want to add driver mutations
# Add statistical tests
# Add other PRSs

# Required data
# m6A peaks
# Somatic mutations [P value]
# - PGA
# - Driver Mutations [Optional]
# - Number of SNVs [Optional]
# - Clonality [Potential]
# Germline [P value]
# - PRS(s)
# Clinical covariate [P value]
# - ISUP
# - Age
# - BCR
# - IDC
# - Metastasis
# Hypoxia [P value]
# Complementary datasets
# - RNA
# - Protein
# - H3K27ac
# - Methylation
# - Mutation
# - Germline European


# Add Missing Samples Function --------------------------------------------

add_missing_samples <- function(df){
  missing_samples <- setdiff(all_samples, df$Sample)
  missing_samples_mat <- matrix(NA, nrow = length(missing_samples), ncol = ncol(df))
  missing_samples_mat <- data.frame(missing_samples_mat)
  colnames(missing_samples_mat) <- colnames(df)
  missing_samples_mat$Sample = missing_samples
  df = rbind(df, missing_samples_mat)
  return(df)
}

# Loading Data ------------------------------------------------------------

# All Samples
all_samples <- all_samples[!all_samples %in% exclude_samples]

# Number of Peaks
metpeak_peaks <- do.call(rbind, lapply(all_samples, function(this_sample){
  peaks <- read.delim(
    paste0("F_MeTPeak_Peaks_SS/", this_sample, "/peak.bed")
  )
  data.frame(
    "Sample" = this_sample,
    "NPeaks" = nrow(peaks)
  )
}))

# Somatic
driver_mutations <- read.delim("Summary_Tables/Bhandari_2019_Mutations/20220421-165034_drivers_m6A_enzyme_mutations.tsv", header = T)
mutational_features <- read.delim("/cluster/home/helenzhu/code/snakemake_M6A_40_PaperFigures/Bhandari_Supplementary/STable4_CPCGENE.Data.txt", header = T)
mutational_features <- mutational_features[
  mutational_features$patient_id %in% all_samples,
  c("patient_id", "PGA", "snv_count", "gr_count")
]
colnames(mutational_features) <- c("Sample", "PGA", "SNVs", "GRs")
mutational_features <- add_missing_samples(mutational_features)

# Germline
# Schumacher PRS
# schumacher <- read.csv("~/Cluster_Helen/Snakemake_Germline_Plink/PRS/Houlahan_PRS_STable1.csv", header = T)
# schumacher <- schumacher[,c("sample", "PRS")]
# schumacher <- schumacher[schumacher$sample %in% all.euro.samples,]
# colnames(schumacher)[colnames(schumacher) == "sample"] <- "Sample"
# schumacher <- add_missing_samples(schumacher)

# Conti PRS
all_PRS <- readRDS("~/Cluster_Helen/Snakemake_Germline_Plink/PRS_Nicole/2023-06-22_Houlahan_CPCG_External_UK_genotype-01_1KG-imputed_processed_prs_data.Rds")
conti <- all_PRS[["conti.multiethnic.prs"]][["per.sample.prs.data"]]
conti$Sample <- gsub("[_].*", "", rownames(conti))
conti <- conti[,c( "Sample", "prs.weighted.sum")]
colnames(conti)[colnames(conti) == "prs.weighted.sum"] <- "PRS"
conti <- conti[conti$Sample %in% all_samples,]
conti <- add_missing_samples(conti)

# Hypoxia
# hypoxia <- read.delim("Summary_Tables/Hypoxia/hypoxia_scores_buffa.txt", header = T)
# hypoxia <- hypoxia[hypoxia$patient.IDs %in% all_samples, c("patient.IDs", "score")]
# colnames(hypoxia)[colnames(hypoxia) == "patient.IDs"] <- "Sample"
# hypoxia <- add_missing_samples(hypoxia)

# Clinical covariates (Age, T Categ, PSA, ISUP/Gleason, BCR)
clinical <- read.delim("Summary_Tables/2019-08-20_500pg_SupTable1_feature_by_patient.tsv", header = T)
clinical$age <- ifelse(is.na(clinical$age_at_diagnosis), clinical$age_at_treatment, clinical$age_at_diagnosis)
clinical <- clinical[
  clinical$sample %in% all_samples,
  c("sample", "age", "pre_treatment_psa", "summary_t", "pathologic_isup_grade")
]
colnames(clinical) <- c("Sample", "Age", "PSA", "TCategory", "ISUP")
clinical <- add_missing_samples(clinical)

# IDC
idc <- read.delim("Summary_Tables/2021-05-07_m6A_idc_or_cribiform.txt", header = T)
clinical$IDC <- as.numeric(idc$idc_or_cribiform[match(clinical$Sample, idc$patient_id)])

# BCR & Metastasis
# TODO: Change colnames
print(load("SummaryFilesNew/bcr_data.rsav"))
bcr_data$sample <- gsub("[-].*", "", rownames(bcr_data))
bcr_data <- bcr_data[bcr_data$sample %in% all_samples,]
bcr_data <- unique(bcr_data)
colnames(bcr_data)[colnames(bcr_data) == "sample"] <- "Sample"

clinical$BCR <- as.numeric(bcr_data$bcr[match(clinical$Sample, bcr_data$Sample)])
clinical$Metastasis <- as.numeric(bcr_data$mets[match(clinical$Sample, bcr_data$Sample)])

# Complementary datasets
complementary_dataset <- data.frame("Sample" = all_samples)
methylation <- readLines("Summary_Tables/methylation_samples.txt")
h3k27ac <- readLines("Summary_Tables/H3K27Ac_samples.txt")
germline_european <- readLines("Summary_Tables/germline_samples.txt")
bulk_rnaseq <- read.csv("Summary_Tables/Chen2019_clinical.csv", header = T)
bulk_rnaseq <- grep("^CPCG", bulk_rnaseq[,"ID"], value = T)
protein <- read.csv("Summary_Tables/Sinha2019_clinical.csv", header = T, skip = 2)[,"Patient_ID"]

complementary_dataset$Methylome <- ifelse(complementary_dataset$Sample %in% methylation, 1, 0)
complementary_dataset$H3K27ac <- ifelse(complementary_dataset$Sample %in% h3k27ac, 1, 0)
complementary_dataset$Germline <- ifelse(complementary_dataset$Sample %in% germline_european, 1, 0)
complementary_dataset$RNAseq <- ifelse(complementary_dataset$Sample %in% bulk_rnaseq, 1, 0)
complementary_dataset$Protein <- ifelse(complementary_dataset$Sample %in% protein, 1, 0)
complementary_dataset$CNA <- ifelse(complementary_dataset$Sample %in% rownames(driver_mutations), 1, 0)
complementary_dataset$SSM <- ifelse(complementary_dataset$Sample %in% rownames(driver_mutations)[is.na(driver_mutations$SSM.SPOP)], 0, 1)
complementary_dataset <- complementary_dataset[,c("Sample", "Germline", "CNA", "SSM", "Methylome", "H3K27ac", "RNAseq", "Protein")]

# Clinical ----------------------------------------------------------------

clinical_colour_scheme = c(
  force.colour.scheme(x = c("<40", "40 - 50", "50 - 65", "65 - 70", ">= 70"), scheme = "age.categorical.prostate"),
  force.colour.scheme(x = c("t1", "t2", "t3", "t4"), scheme = "clinicalt3"),
  force.colour.scheme(x = c("0 - 9.9", "10 - 19.9", ">= 20"), scheme = "psa.categorical"),
  force.colour.scheme(x = c("1", "2", "3", "4", "5"), scheme = "isup.grade"),
  c("white", "darkred"),
  c("white", "blue4"),
  c("white", "deeppink3")
)

# Encoding clinical data
clinical_encoded <- clinical
clinical_encoded$Age <- cut(clinical_encoded$Age, breaks = c(0, 40, 50, 65, 70, 100))
levels(clinical_encoded$Age) <- c("<40", "40 - 50", "50 - 65", "65 - 70", ">= 70")
clinical_encoded$TCategory <- factor(
  substr(clinical_encoded$TCategory, 1, 2),
  levels = c("T1", "T2", "T3", "T4")
)
clinical_encoded$PSA <-  cut(clinical_encoded$PSA, breaks = c(0, 10, 20, 30))
levels(clinical_encoded$PSA) <- c("0 - 9.9", "10 - 19.9", ">= 20")

# Order (Age, TCategory, PSA, ISUP)
# Age: 1 - 5
# T Category: 6 - 9
# PSA: 10 - 12
# ISUP: 13 - 17
# IDC: 18 - 19
# BCR: 20 - 21
# Met: 22 - 23

clinical_plotting <- clinical_encoded[
  c("ISUP", "TCategory", "PSA", "Age", "IDC", "BCR", "Metastasis")
]
rownames(clinical_plotting) <- clinical_encoded$Sample

clinical_plotting$Age <- as.numeric(clinical_plotting$Age)
clinical_plotting$TCategory <- as.numeric(clinical_plotting$TCategory) + 5
clinical_plotting$PSA <- as.numeric(clinical_plotting$PSA) + 9
clinical_plotting$ISUP <- clinical_plotting$ISUP + 12
clinical_plotting$IDC <- clinical_plotting$IDC + 18
clinical_plotting$BCR <- clinical_plotting$BCR + 20
clinical_plotting$Metastasis <- clinical_plotting$Metastasis + 22

# Plotting ----------------------------------------------------------------

# Number of peaks barplot
metpeak_peaks <- metpeak_peaks[order(metpeak_peaks$NPeaks, decreasing = T),]
metpeak_peaks$Sample <- factor(metpeak_peaks$Sample, levels = metpeak_peaks$Sample)
peak_barplot <- create.barplot(
  NPeaks/1000 ~ Sample,
  metpeak_peaks,
  # Barplot formatting
  col = "transparent",
  # Y axis formatting
  ylimits = c(0, 12),
  yat = seq(0, 12, 6),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("Peaks (10"^3*")")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

# Number of peaks lineplot
peak_lineplot <- create.scatterplot(
  NPeaks/1000 ~ Sample,
  metpeak_peaks,
  # Scatterplot formatting
  type = "h",
  # Y axis formatting
  ylimits = c(0, 12),
  yat = seq(0, 12, 6),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("Peaks (10"^3*")")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

# Sets the order for the number of samples
sample_order <- levels(metpeak_peaks$Sample)

# Mutational features barplot
mutational_features$Sample <- factor(mutational_features$Sample, levels = sample_order)

sv_barplot <- create.barplot(
  GRs/1000 ~ Sample,
  mutational_features,
  # Barplot formatting
  col = "transparent",
  # Y axis formatting
  ylimits = c(0, 2),
  yat = seq(0, 2, 1),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("GRs (10"^3*")")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

sv_lineplot <- create.scatterplot(
  GRs/1000 ~ Sample,
  mutational_features,
  # scatterplot formatting
  type = "h",
  # Y axis formatting
  ylimits = c(0, 2),
  yat = seq(0, 2, 1),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("GRs (10"^3*")")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

SNVs_barplot <- create.barplot(
  SNVs/1000 ~ Sample,
  mutational_features,
  # Barplot formatting
  col = "transparent",
  # Y axis formatting
  ylimits = c(0, 6),
  yat = seq(0, 6, 2),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("SNVs (10"^3*")")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

snv_lineplot <- create.scatterplot(
  SNVs/1000 ~ Sample,
  mutational_features,
  # scatterplot formatting
  type = "h",
  # Y axis formatting
  ylimits = c(0, 6),
  yat = seq(0, 6, 2),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("SNVs (10"^3*")")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

PGA_barplot <- create.barplot(
  PGA ~ Sample,
  mutational_features,
  # Barplot formatting
  col = "transparent",
  # Y axis formatting
  ylimits = c(0, 40),
  yat = seq(0, 40, 20),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("PGA")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

PGA_lineplot <- create.scatterplot(
  PGA ~ Sample,
  mutational_features,
  # scatterplot formatting
  type = "h",
  # Y axis formatting
  ylimits = c(0, 40),
  yat = seq(0, 40, 20),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("PGA")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

# # Hypoxia plot
# hypoxia$Sample <- factor(hypoxia$Sample, levels = sample_order)
# hypoxia_lollipop_plot <- create.scatterplot(
#   score ~ Sample,
#   hypoxia,
#   # Scatterplot formatting
#   type = c("p", "h"),
#   abline.h = c(0),
#   cex = 0.5,
#   # Y axis formatting
#   ylimits = c(-35, 35),
#   yat = c(-30, 0, 30),
#   # General formatting
#   xlab.cex = 0,
#   ylab.cex = 0.8,
#   ylab.label = expression(bold("Hypoxia")),
#   xaxis.cex = 0,
#   yaxis.cex = 0.8,
#   xaxis.rot = 90,
#   xaxis.tck = 0,
#   yaxis.tck = 0
# )

# PRS plot
conti$Sample <- factor(conti$Sample, levels = sample_order)
prs_lollipop_plot <- create.scatterplot(
  PRS ~ Sample,
  conti,
  # Scatterplot formatting
  type = c("p", "h"),
  cex = 0.5,
  # Y axis formatting
  ylimits = c(20, 25),
  yat = c(20, 25),
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0.8,
  ylab.label = expression(bold("PRS")),
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0
)

# Clinical covariate plot
clinical_heatmap <- create.heatmap(
  t(clinical_plotting[all_samples,]),
  # Heatmap parameters
  clustering.method = "none",
  same.as.matrix = T,
  # Heatmap coloring
  colour.scheme = clinical_colour_scheme,
  at = seq(0, length(clinical_colour_scheme)) + 0.5,
  total.colours = length(clinical_colour_scheme) + 1,
  print.colour.key = F,
  axes.lwd = 1,
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0,
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0,
  # Other
  yaxis.lab = c("ISUP", "T Category", "PSA", "Age", "IDC/CA", "BCR", "Metastasis")
)

# Clinical barplot
clinical_barplot_data <- do.call(
  rbind.data.frame,
  lapply(c("ISUP", "TCategory", "PSA", "Age", "IDC", "BCR", "Metastasis"), function(ID) {
    res <- data.frame(table(clinical_encoded[,ID]))
    res$ID <- ID
    return(res)
  }
))
clinical_barplot_data$Var1 <- as.character(clinical_barplot_data$Var1)
clinical_barplot_data <- clinical_barplot_data[clinical_barplot_data$Var1 != "0",]
clinical_barplot_data$Var1[clinical_barplot_data$ID %in% c("IDC", "BCR", "Metastasis")] <- clinical_barplot_data$ID[clinical_barplot_data$ID %in% c("IDC", "BCR", "Metastasis")]
clinical_barplot_data$ID <- factor(
  clinical_barplot_data$ID,
  levels = rev(c("ISUP", "TCategory", "PSA", "Age", "IDC", "BCR", "Metastasis"))
)
clinical_barplot_data$Var1 <- factor(clinical_barplot_data$Var1, levels = clinical_barplot_data$Var1)
clinical_barplot_data_colour_scheme <- c(
  force.colour.scheme(x = c("1", "2", "3", "4", "5"), scheme = "isup.grade"),
  force.colour.scheme(x = c("t1", "t2", "t3", "t4"), scheme = "clinicalt3"),
  force.colour.scheme(x = c("0 - 9.9", "10 - 19.9", ">= 20"), scheme = "psa.categorical"),
  force.colour.scheme(x = c("<40", "40 - 50", "50 - 65", "65 - 70", ">= 70"), scheme = "age.categorical.prostate"),
  c("darkred", "blue4", "deeppink3")
)
clinical_barplot <- create.barplot(
  ID ~ Freq,
  clinical_barplot_data,
  groups = clinical_barplot_data$Var1,
  col = clinical_barplot_data_colour_scheme,
  stack = TRUE,
  plot.horizontal = T,
  # Formatting
  xaxis.cex = 0,
  yaxis.cex = 0,
  xlab.cex = 0,
  ylab.cex = 0,
  xlimits = c(0, 148),
  xaxis.tck = 0,
  yaxis.tck = 0
)

# Additional datasets plot
complementary_dataset_plotting <- complementary_dataset[!colnames(complementary_dataset) == "Sample"]
rownames(complementary_dataset_plotting) <- complementary_dataset$Sample
complementary_datasets_heatmap <- create.heatmap(
  t(complementary_dataset_plotting[all_samples,]),
  # Heatmap parameters
  clustering.method = "none",
  same.as.matrix = T,
  # Heatmap coloring
  colour.scheme = c("white", "black"),
  at = c(-0.5, 0.5, 1.5),
  total.colours = 3,
  print.colour.key = F,
  axes.lwd = 1,
  # General formatting
  xlab.cex = 0,
  ylab.cex = 0,
  xaxis.cex = 0,
  yaxis.cex = 0.8,
  xaxis.rot = 90,
  xaxis.tck = 0,
  yaxis.tck = 0,
  # Other
  yaxis.lab = colnames(complementary_dataset_plotting)
)

# Additional datasets barplot
complementary_dataset_barplot_data <- data.frame(
  "Type" = colnames(complementary_dataset_plotting),
  "Count" = colSums(complementary_dataset_plotting)
)
complementary_dataset_barplot_data$Type = factor(
  complementary_dataset_barplot_data$Type,
  levels = rev(complementary_dataset_barplot_data$Type)
)
complementary_datasets_barplot <- create.barplot(
  Type ~ Count,
  complementary_dataset_barplot_data,
  # Barplot
  plot.horizontal = TRUE,
  # Formatting
  xaxis.cex = 1,
  yaxis.cex = 0,
  xlab.cex = 0,
  ylab.cex = 0,
  xlimits = c(0, 148),
  xat = c(0, 148),
  xaxis.tck = 0,
  yaxis.tck = 0
)

# Legend
covariate.legend <- list(
  legend = list(
    colours = c(force.colour.scheme(x = c("1", "2", "3", "4", "5"), scheme = "isup.grade"), "slategray"),
    labels = c("1", "2", "3", "4", "5", "N/A"),
    title = expression(bold(underline('ISUP Grade'))),
    lwd = 0.5
  ),
  legend = list(
    colours = c(force.colour.scheme(x = c("0 - 9.9", "10 - 19.9", ">= 20"), scheme = "psa.categorical"), "slategray"),
    labels = c("0 - 10", expression("" > "10 - 20"), expression("" > "20"), "N/A"),
    title = expression(bold(underline('PSA (ng/mL)'))),
    lwd = 0.5
  ),
  legend = list(
    colours = force.colour.scheme(x = c("t1", "t2", "t3", "t4"), scheme = "clinicalt3"),
    labels = c("T1", "T2", "T3", "T4"),
    title = expression(bold(underline('T Category'))),
    lwd = 0.5
  ),
  legend = list(
    colours = c(force.colour.scheme(x = c("<40", "40 - 50", "50 - 65", "65 - 70", ">= 70"), scheme = "age.categorical.prostate"), "slategray"),
    labels = c(expression("" < "40"), "40 - 50", "50 - 65", "65 - 70", expression("" >= "70"), "N/A"),
    title = expression(bold(underline('Age'))),
    lwd = 0.5
  ),
  legend = list(
    colours = c("darkred", "slategray"),
    labels = c("TRUE", "N/A"),
    title = expression(bold(underline('IDC/CA'))),
    lwd = 0.5
  ),
  legend = list(
    colours = c("blue4", "slategray"),
    labels = c("TRUE", "N/A"),
    title = expression(bold(underline('BCR'))),
    lwd = 0.5
  ),
  legend = list(
    colours = c("deeppink3", "slategray"),
    labels = c("TRUE", "N/A"),
    title = expression(bold(underline('Metastasis'))),
    lwd = 0.5
  )
)

side.legend <-legend.grob(
  legends = covariate.legend,
  label.cex = 0.7,
  title.cex = 0.7,
  title.just = 'left',
  title.fontface = 'bold',
  size = 2,
  layout = c(2, 4),
  x = -0.85,
  y = 0.65
)

# Creating the multiplanelplot
plot_objects <- list(
  # peak_barplot,
  peak_lineplot,
  # sv_barplot,
  sv_lineplot,
  # SNVs_barplot,
  snv_lineplot,
  # PGA_barplot,
  PGA_lineplot,
  # hypoxia_lollipop_plot,
  prs_lollipop_plot,
  clinical_heatmap,
  clinical_barplot,
  complementary_datasets_heatmap,
  complementary_datasets_barplot
)

compiled_plot <- create.multipanelplot(
  plot.objects = plot_objects,
  # Layout
  layout.width = 2,
  layout.height = length(plot_objects) - 2,
  layout.skip = c(
    F, T,
    F, T,
    F, T,
    F, T,
    F, T,
    F, F,
    F, F
  ),
  # Spacing
  x.spacing = 0.5,
  y.spacing = c(-3, -3, -3, -3, -3, -3, -3),
  # Object dimensions
  plot.objects.heights = c(1, 1, 1, 1, 1, 1.3, 1.3),
  plot.objects.widths = c(1, 0.15),
  # Axes padding
  ylab.axis.padding = -6,
  # General padding
  left.padding = 0,
  right.padding = 0,
  top.padding = 0,
  bottom.padding = 0,
  # Legend
  legend = list(
    right = list(
      x = 0.8,
      y = 1,
      fun = side.legend
    )
  )
)


# Saving the plot
filename = '~/figures/091_OverviewFigure.pdf'
pdf(filename, width = 14, height = 7.5)
print(compiled_plot)
dev.off()


# Statistical Tests -------------------------------------------------------

# Combining data
combined_data <- Reduce(function(df1, df2) merge(df1, df2, by = "Sample", all = T),
                        list(metpeak_peaks, mutational_features, conti, clinical)) #  hypoxia,

# Correlation results
todo <- c("PGA", "SNVs", "GRs", "PRS","Age", "PSA") #  "score",
correlation_results <- do.call(rbind.data.frame, lapply(todo, function(id){
  res = cor.test(combined_data[,id], combined_data[,"NPeaks"], method = "spearman")
  data.frame(
    "ID" = id,
    "EffectSize" = res$estimate,
    "P.value" = res$p.value,
    "Test" = "spearman"
  )
}))

# T-test results
idc_res <- t.test(combined_data$NPeaks[combined_data$IDC == 1], combined_data$NPeaks[combined_data$IDC == 0])
ttest_results <- data.frame(
  "ID" = "IDC",
  "EffectSize" = mean(combined_data$NPeaks[combined_data$IDC == 1], na.rm = T)/mean(combined_data$NPeaks[combined_data$IDC == 0], na.rm = T),
  "P.value" = idc_res$p.value,
  "Test" = "t-test"

)

# Anova
todo <- c("ISUP", "TCategory")
aov_results <- do.call(rbind.data.frame, lapply(todo, function(id){
  aov_res <- aov(NPeaks ~ combined_data[,id], combined_data)
  data.frame(
    "ID" = id,
    "EffectSize" = NA,
    "P.value" = summary(aov_res)[[1]][["Pr(>F)"]][1],
    "Test" = "ANOVA"
  )
}))

# Survival Analysis
# Tried it with median dichotomization, didn't work
bcr_combined <- merge(metpeak_peaks, bcr_data, by = "Sample")
bcr_combined$dichotomy <- ifelse(bcr_combined$NPeaks > median(bcr_combined$NPeaks), 1, 0)
bcr_cox <- coxph(Surv(time_to_bcr, bcr) ~ NPeaks, data = bcr_combined)
met_cox <- coxph(Surv(time_to_mets, mets) ~ NPeaks, data = bcr_combined)

surv_results <- data.frame(
  "ID" = c("BCR", "Mets"),
  "EffectSize" = c(summary(bcr_cox)$conf.int[1,1], summary(met_cox)$conf.int[1,1]),
  "P.value" = c(summary(bcr_cox)$waldtest['pvalue'], summary(met_cox)$waldtest['pvalue']),
  "Test" = "CoxPH"
)

all_results <- rbind(correlation_results, ttest_results, aov_results, surv_results)
clinical_results <- all_results[all_results$ID %in% c("PGA", "Age", "PSA", "IDC", "ISUP", "TCategory", "BCR", "Mets"),] # "PRS",

# Adding an association dotmap --------------------------------------------

key.func = function(x) { 0.1 + abs(x) ; }
key.sizes = c(-2, -1, -0.5, 0, 0.5, 1, 2)

key.tmp = list(
  space = 'right',
  points = list(
    cex =  key.func(key.sizes),
    col = c(rep("dodgerblue2", 3), "transparent", rep("darkorange1", 3)),
    pch = 19
  ),
  text = list(
    lab = as.character(key.sizes),
    cex = 1,
    adj = 1
  ),
  padding.text = 3,
  background = 'white'
)

clinical_results <- clinical_results[
  match(c("PGA", "Age", "PSA", "ISUP", "TCategory", "IDC", "BCR", "Mets"), clinical_results$ID), # "PRS",
]
clinical_results$lab <-  c(
  # expression("PRS "*rho),
  expression("PGA "*rho),
  expression("Age "*rho),
  expression("PSA "*rho),
  expression("ISUP"),
  expression("T Category"),
  expression("IDC FC"),
  expression("BCR HR"),
  expression("Metastasis HR")
)

# Clinical dotmap
clinical_dotmap <- create.dotmap(
  x = matrix(clinical_results$EffectSize),
  bg.data = matrix(-log10(clinical_results$P.value)),
  # Background formatting
  colour.scheme = c('white','black'),
  bg.alpha = 1,
  at = c(0, 1, 2),
  colourkey = TRUE,
  colourkey.labels.at = c(0.5, 1.5),
  colourkey.labels = c(
    "N.S.",
    "P < 0.1"
  ),
  # Spot size
  spot.size.function = key.func,
  # lwd
  col.lwd = 1,
  row.lwd = 1,
  lwd = 1,
  # Formatting
  xaxis.cex = 0,
  yaxis.cex = 1,
  yaxis.lab = clinical_results$lab,
  xlab.cex = 0,
  ylab.cex = 0,
  xaxis.tck = 0,
  yaxis.tck = 0,
  # Key
  key = key.tmp,
  key.top = 1,
  # NA
  na.spot.size = 3,
  na.pch = 4,
  na.spot.size.colour = 'black',
)

# Printing plot
filename <- "~/figures/091_overview_dotmap.pdf"
pdf(filename, width = 3.5, height = 4)
print(clinical_dotmap)
dev.off()

# Boxplots and scatterplots -----------------------------------------------

plotting_data <- combined_data
plotting_data$IDC <- ifelse(plotting_data$IDC == 1, "Yes", "No")

pval_peaks_idc <- formatC(clinical_results$P.value[clinical_results$ID == "IDC"], 2)
delta_peaks_idc <- mean(combined_data$NPeaks[combined_data$IDC == 1], na.rm = T) - mean(combined_data$NPeaks[combined_data$IDC == 0], na.rm = T)
delta_peaks_idc <- formatC(delta_peaks_idc, 3)

idc_boxplot <- create.boxplot(
  NPeaks ~ IDC,
  plotting_data,
  # Boxplot formatting
  add.stripplot = T,
  # Formatting
  xaxis.cex = 1,
  xaxis.lab = c("No IDC/CA", "IDC/CA"),
  yaxis.cex = 1,
  yat = seq(0, 12000, 4000),
  ylimits = c(0, 12000),
  ylab.cex = 1,
  ylab.label = "Peaks",
  xlab.cex = 0,
  # Adding P value
  add.text = TRUE,
  text.labels = as.expression(
    bquote(
      atop(
      "P < " ~ .(pval_peaks_idc),
      Delta ~ "Peaks =" ~ .(delta_peaks_idc)
  ))),
  text.x = 2,
  text.y = 1100,
  text.anchor = 'centre',
  text.col = 'black',
  text.cex = 0.8,
  text.fontface = 1,
)

prs_scatterplot <- create.scatterplot(
  PRS ~ NPeaks,
  plotting_data,
  # Formatting
  xaxis.cex = 1,
  yaxis.cex = 1,
  xaxis.tck = 0,
  yaxis.tck = 0,
  xlimits = c(0, 12000),
  xat = seq(0, 12000, 4000),
  ylimits = c(20, 25),
  yat = seq(20, 25, 1),
  ylab.cex = 1,
  xlab.cex = 1,
  ylab.label = "Conti PRS",
  xlab.label = "Peaks",
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = get.corr.key(
          # two vectors of the same length
          x = plotting_data$PRS,
          y = plotting_data$NPeaks,
          # specify what to include in the key
          label.items = c('spearman','spearman.p'),
          # format key
          alpha.background = 0,
          key.cex = 1
        )
      ),
      # place key
      x = 0.6,
      y = 0.2,
      corner = c(0,1)
    )
  )
)

# Printing plot
filename <- "~/figures/091_oveview_IDC_boxplot.pdf"
pdf(filename, width = 4, height = 3)
print(idc_boxplot)
dev.off()

filename <- "~/figures/091_oveview_PRS_scatterplot.pdf"
pdf(filename, width = 4, height = 4)
print(prs_scatterplot)
dev.off()
