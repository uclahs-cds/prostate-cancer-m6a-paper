source("/cluster/home/helenzhu/code/snakemake_M6A_36_m6AQTLDB/Rscripts/000_HEADER.R")

library(vcfR)
library(GenomicRanges)

# Preamble ----------------------------------------------------------------
# This looks at the allele frequencies of the variants

# Loading Data ------------------------------------------------------------

# Variants
vcf_file = "Database.m6A.analysis/Database.Sites.vcf.gz"
variants = read.vcfR(vcf_file)
info = data.frame(getFIX(variants))
info$POS = strtoi(info$POS)

# Allele Frequency
maf_file = "Database.m6A.analysis/Database.Sites.frq"
maf = read.table(
  maf_file, 
  header = F, 
  sep = "\t",
  skip = 1,
  stringsAsFactors = F)
colnames(maf) = c("CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE:FREQ1", "ALLELE:FREQ2")

# Peaks
peak.file = "/cluster/home/helenzhu/Cluster_Helen/Snakemake_M6A_hg38/H_PeakCounts/MeTPeak.V9_MetricVoting_OptimizationFunction_MaxGapOnly.bed6"
peaks = read.delim(peak.file, header = F)
colnames(peaks) = c("chr", "start", "end", "name", "score", "strand")

# Analysis ----------------------------------------------------------------

# Check
all(info$REF %in% c("A", "T") | info$ALT %in% c("A", "T"))

# Calculating allele frequency
all.results = merge(info, maf, by.x = c("CHROM", "POS"), by.y = c("CHROM", "POS"))
all.results$A.Freq = ifelse(all.results$REF %in% c("A", "T"), all.results[,"ALLELE:FREQ1"], all.results[,"ALLELE:FREQ2"])

# Finding Overlap with Peaks
peaks.gr = GRanges(
  seqnames = peaks$chr,
  IRanges(
    start = peaks$start + 1, # Bed files are base 0
    end = peaks$end
  )
)

variants.gr = GRanges(
  seqnames = all.results$CHROM,
  IRanges(
    start = all.results$POS,
    end = all.results$POS
  )
)
ovl = findOverlaps(peaks.gr, variants.gr)

# Peak variants
peak.variants = cbind(all.results[subjectHits(ovl),], peaks[queryHits(ovl),])
peak.variants$gene = gene_conversion[gsub("[:].*", "", peak.variants$name)]
peak.variants = peak.variants[(peak.variants$REF == "A" & peak.variants$strand == "+") |
                                (peak.variants$REF == "T" & peak.variants$strand == "-"),]

# Annotating all results
all.results$PeakStatus = ifelse(all.results$ID %in% peak.variants$ID, 1, 0)

# Plotting ----------------------------------------------------------------

barplot.data = data.frame(
  "AF" = cut(all.results$A.Freq, breaks = seq(0, 1, 0.1), include.lowest = T),
  "PeakStatus" = all.results$PeakStatus
)
barplot.data = data.frame(
  table(barplot.data$AF, barplot.data$PeakStatus)
)
colnames(barplot.data) = c("AF", "PeakStatus", "Count")
barplot.data$Sum = ave(barplot.data$Count, barplot.data$PeakStatus, FUN = sum)
barplot.data$Percent = barplot.data$Count/barplot.data$Sum*100

filename = "~/figures/110_AlleleFrequency.m6A.variants.pdf"
pdf(filename, width = 3.5, height = 3)

create.barplot(
  Percent ~ AF,
  barplot.data,
  groups = barplot.data$PeakStatus,
  col = c("transparent", "red"),
  xaxis.cex = 1,
  yaxis.cex = 1,
  ylab.label = "% SNPs",
  xlab.label = "Allele Frequency",
  ylab.cex = 1,
  xlab.cex = 1,
  xat = seq(0, 10, 2) + 0.5,
  xaxis.lab = seq(0, 1, 0.2),
  xaxis.tck = 0,
  yaxis.tck = 0,
  box.ratio = 10,
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = list(
          points = list(
            col = 'black',
            pch = c(22, 22),
            cex = 1.5,
            # reverse order to match stacked bar order
            fill = c('transparent', 'red')
          ),
          text = list(
            # reverse order to match stacked bar order
            lab = c(expression("All m"^"6"*"A variants"),
                    expression("PCa m"^"6"*"A variants")),
            cex = 0.7
          ),
          padding.text = 0.7,
          cex = 0.5
        )
      ),
      x = 0.0,
      y = 0.98
    )
  )
)

dev.off()



# Writing Output Tables ---------------------------------------------------

filename = "Database.m6A.analysis/All.m6A.Site.Variants.tsv"
write.table(
  all.results,
  file = filename,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

filename = "Database.m6A.analysis/m6A.Peak.Variants.tsv"
write.table(
  peak.variants,
  file = filename,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
system(paste0("cp ", filename, " ~"))