# Associated code for "The Landscape of N6-Methyladenosine in Localized Primary Prostate Cancer"
## Notes shown below where applicable
## Pipelines

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

**F_DataNormalization**
* Normalizing and adjusting IP data by Input data

**G_AllelicImbalance**
* Generating allelic-specific reads for Figure 3A

**H_m6AQTLAnalysis**
* Quality control of germline variants
* Testing Plink/PEER covariates
* Generating input files and running MatrixQTL

**I_HypoxiaCellLines**
* Processing hypoxia cell line meRIP-seq data

**J_TumourNormal**
* Processing tumour normal meRIP-seq data

**SAC-seqProcessing**
* Data processed following code from https://github.com/shunliubio/m6A-SAC-seq_new_protocol

## Figures

**Figure 3**
* The code for Figure 3B is provided in the code for Supplementary Figure 4D-F
* The code for Figure 3C is provided in the code for Figure 1A
* The code for Figure 3E-F is provided in the code for SupplementaryFigure 4G

**Figure 4**
* Code provided for Figure 4A

**Figure 5**
* Code provided for Figure 5B, E, J, K

## Supplementary Figures

**Supplementary Figure 1**
* Supplementary Figures 1A-B were drawn rather than coded
* Supplementary Figure 1E - Data generated from D_RNAseqGermlineVariants pipeline

**Supplementary Figure 2**
* Supplementary Figure 2C-D - Data generated from C_SubsamplingPeakAnalysis pipeline
* Supplementary Figures 2G-J are addressed by Figure1/B-E_GeneProfilePlots.R

**Supplementary Figure 4**
* Code provided for bigWig processing for producing Figure 4 D & E

**Supplementary Figure 5**
* Supplementary Figures 5A & 5B are from the ConsensusClusterPlus package

**Supplementary Figure 6**
* code for Supplementary Figure 6C is in Figure1/A_Overview.R (see idc_boxplot)
* code for Supplementary Figure 6D is in Figure1/A_Overview.R (see clinical_dotmap)

**Supplementary Figure 7**
* Supplementary Figures 7A & 7J were drawn rather than coded
* The code for Supplementary Figure 7I is provided in the code for Figure 1A
