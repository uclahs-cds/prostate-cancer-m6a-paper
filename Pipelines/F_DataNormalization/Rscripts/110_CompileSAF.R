


# Preamble ----------------------------------------------------------------
# This creates SAF files from BeD6 files


# Loading Data ------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
input.file = args[1]
output.file = args[2]


# Output ------------------------------------------------------------------

# Reading Table
bed6 = read.table(input.file, header = F, sep = "\t", stringsAsFactors = F)

# Column Names
colnames(bed6) = c("Chr", "Start", "End", "GeneID", "Score", "Strand")

# Converting base 0 to base 1
bed6$Start = bed6$Start+1

# SAF format
saf = bed6[,c("GeneID", "Chr", "Start", "End", "Strand")]

# No Scientific Notation in Output
options(scipen = 999)

# Output
write.table(
  saf,
  file = output.file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)