# Code for the prostate cancer m6A profiling study

## Pipelines
Code for processing data
TODO: Pipeline for QC compilation
TODO: Pipeline for m6A-SAC-seq data

**A_ProcessingRawData**
* Read alignment and quality control
* Generation of BigWig files

**B_PeakCalling**
* Calling peaks using MeTPeak and exomePeak
* Calling motifs using HOMER

**C_SubsamplingPeakAnalysis**
* Subsampling reads from the IP and Input library and evaluating their impact on peak calling
* Results used to generate Supplementary Figure 2C-D

**D_RNAseqGermlineVariants**
* Calling germline variants using the Input library and comparing them to the variants called in Houlahan *et al.,* 2019
* Results used to generate Supplementary Figure 1E

**E_HistogramZoo**
* Calling joint peaks using HistogramZoo
* Normalizing and adjusting IP data by Input data

**F_AllelicImbalance**
* Generating allelic-specific reads for Figure 3A

**G_m6AQTLAnalysis**

**H_HypoxiaCellLines**
* Processing hypoxia cell line meRIP-seq data

**I_TumourNormal**
* Processing tumour normal meRIP-seq data


## Figures

**Figure 1**
* TODO: Add the code for Figure 1G

**Figure 3**
* The code for Figure 3B is provided in the code for Supplementary Figure 4D-F
* The code for Figure 3C is provided in the code for Figure 1A
* The code for Figure 3E-F is provided in the code for SupplementaryFigure 4G

## Supplementary Figures

**Supplementary Figure 1**
* Supplementary Figures 1A-B were drawn rather than coded
* TODO: Supplementary Figure 1C: QC data not provided, compiled in other scripts
* Supplementary Figure 1E - Data generated from D_RNAseqGermlineVariants pipeline

**Supplementary Figure 2**
* Supplementary Figure 2C-D - Data generated from C_SubsamplingPeakAnalysis pipeline
* Supplementary Figures 2G-J are addressed by Figure1/B-E_GeneProfilePlots.R

**Supplementary Figure 4**
* Supplementary Figures 4A & 4J were drawn rather than coded
* The code for Supplementary Figure 4I is provided in the code for Figure 1A

**Supplementary Figure 5**
