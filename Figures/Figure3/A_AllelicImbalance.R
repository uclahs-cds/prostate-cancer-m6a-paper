source("/cluster/home/helenzhu/code/snakemake_M6A_36_m6AQTLDB/Rscripts/000_HEADER.R")



# Preamble ----------------------------------------------------------------
# This tests for allelic imbalance in the 36 variants
# Both Input & IP

# Loading Data ------------------------------------------------------------

# Loading Peak Variants
peak.variants.file = "Database.m6A.analysis/m6A.Peak.Variants.tsv"
peak.variants = read.delim(peak.variants.file, header = T)

# # Loading IP file
# ip.file = "AE_Matrices/ip.counts.tsv"
# ip.counts = read.delim(ip.file, header = T)
# ip.counts = ip.counts[ip.counts$variantID %in% peak.variants$ID,]
# 
# # Loading Input file
# input.file = "AE_Matrices/input.counts.tsv"
# input.counts = read.delim(input.file, header = T)
# input.counts = input.counts[input.counts$variantID %in% peak.variants$ID,]
# 
# # Saving and loading files
# save(ip.counts, input.counts, file = "Database.m6A.analysis/Peak.Variants.AE.rsav")
print(load("Database.m6A.analysis/Peak.Variants.AE.rsav"))
# [1] "ip.counts"    "input.counts"

# Analysis ----------------------------------------------------------------

format.df = function(df){
  df = df[,c("variantID", grep("^CPCG", colnames(df), value = T))]
  df = melt(df, id.vars = "variantID")
  df = df[complete.cases(df),]
  df$Sample = gsub("[.].*", "", as.character(df$variable))
  df$variable = gsub(".*[.]", "", as.character(df$variable))
  df$variantID = as.character(df$variantID)
  df = df[df$variable %in% c("refCount", "altCount"),]
  return(df)  
}

ip.counts = format.df(ip.counts)
input.counts = format.df(input.counts)

# Try running a paired t-test
paired.t = function(snp, df){
  # snp = "rs1043550:chr7:128769171:A:G"
  # df = ip.counts
  df = df[df$variantID == snp,]
  df = dcast(data =  df, formula =  Sample ~ variable, value.var = "value")
  pval = t.test(df$refCount, df$altCount,  paired = T)$p.value
  data.frame(
    "snp" = snp,
    "p.value" = pval,
    "mean.ref" = mean(df$refCount),
    "mean.alt" = mean(df$altCount),
    "foldChange" = mean(df$refCount)/mean(df$altCount)
  )
}

ids = sort(unique(ip.counts$variantID))
ip.results = do.call(rbind, lapply(ids, function(id){ paired.t(snp = id, df = ip.counts)}))
ip.results$q.value = p.adjust(ip.results$p.value, method = "fdr")
colnames(ip.results)[colnames(ip.results) %in% c("p.value", "mean.ref", "mean.alt", "foldChange", "q.value")] = 
  paste0("IP.", colnames(ip.results)[colnames(ip.results) %in% c("p.value", "mean.ref", "mean.alt", "foldChange", "q.value")])

input.results = do.call(rbind, lapply(ids, function(id){ paired.t(snp = id, df = input.counts)}))
input.results$q.value = p.adjust(input.results$p.value, method = "fdr")
colnames(input.results)[colnames(input.results) %in% c("p.value", "mean.ref", "mean.alt", "foldChange", "q.value")] = 
  paste0("Input.", colnames(input.results)[colnames(input.results) %in% c("p.value", "mean.ref", "mean.alt", "foldChange", "q.value")])

all.results = merge(ip.results, input.results, by = "snp")
all.results = merge(peak.variants[,c("ID", "REF", "ALT", "A.Freq", "strand", "name", "gene")], all.results, by.x = "ID", by.y = "snp")

# Making a Dotmap ---------------------------------------------------------

all.results = all.results[order(all.results$IP.q.value, all.results$Input.q.value),]
all.results$tag = paste0(gsub("[:].*", "", all.results$ID), " | ", all.results$gene, " | ", gsub(".*[:]", "", all.results$name))
sig.results = all.results[all.results$IP.q.value < 0.1,]

filename = "~/figures/120_AllelicImbalance.m6A.sites.pdf"
pdf(filename, width = 6, height = 8)

bg.data = -log10(all.results[,c("IP.q.value", "Input.q.value")])
dot.size = log2(all.results[,c("IP.foldChange", "Input.foldChange")])

key.func = function(x) { 0.1 + abs(x) ; }
key.sizes = c(2, -1, -0.5, 0, 0.5, 1, 2)

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
  title = expression("Log"["2"]*"FC"),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)


create.dotmap(
  #
  x = dot.size,
  bg.data = bg.data,
  colour.scheme = c('white','black'),
  filename = NULL,
  # Scaling
  main.cex = 0,
  main = "",
  yaxis.lab = all.results$tag,
  xaxis.lab = c("IP", "Input"),
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  ylab.label = "SNP | Gene | Peak",
  ylab.cex = 1,
  xlab.cex = 0,
  xaxis.tck = 0,
  yaxis.tck = 0,
  # Spot size
  spot.size.function = key.func,
  # lwd
  col.lwd = 1, 
  row.lwd = 1,
  lwd = 1,
  # NA
  na.spot.size = 1, 
  na.pch = 4, 
  na.spot.size.colour = 'black',
  # Key
  key = key.tmp,
  key.top = 1,
  # Background colour
  bg.alpha = 1,
  at = c(0, 1, 2, 8),
  # Background colourkey
  colourkey = TRUE,
  colourkey.labels.at = c(0, 1, 2, 8), 
  colourkey.labels = c(
    expression("1"^" "),
    expression("10"^"-1"),
    expression("10"^"-2"),
    expression("10"^"-8")
  ), 
  description = 'Dotmap created by BoutrosLab.plotting.general',
  resolution = 50
)

dev.off()


filename = "~/figures/120_AllelicImbalance.significant.m6A.sites.pdf"
pdf(filename, width = 6, height = 4)

bg.data = -log10(sig.results[c("IP.q.value", "Input.q.value")])
dot.size = log2(sig.results[c("IP.foldChange", "Input.foldChange")])

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
  title = expression("Log"["2"]*"FC"),
  cex.title = 1,
  padding.text = 3,
  background = 'white'
)


create.dotmap(
  #
  x = dot.size,
  bg.data = bg.data,
  colour.scheme = c('white','black'),
  filename = NULL,
  # Scaling
  main.cex = 0,
  main = "",
  yaxis.lab = sig.results$tag,
  xaxis.lab = c("IP", "Input"),
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  ylab.label = "SNP | Gene | Peak",
  ylab.cex = 1,
  xlab.cex = 1,
  xlab.label = expression("Q value"),
  xaxis.tck = 0,
  yaxis.tck = 0,
  # Spot size
  spot.size.function = key.func,
  # lwd
  col.lwd = 1, 
  row.lwd = 1,
  lwd = 1,
  # NA
  na.spot.size = 1, 
  na.pch = 4, 
  na.spot.size.colour = 'black',
  # Key
  key = key.tmp,
  key.top = 1,
  # Background colour
  bg.alpha = 1,
  at = c(0, 1, 2, 8),
  # Background colourkey
  colourkey = TRUE,
  colourkey.labels.at = c(0, 1, 2, 8), 
  colourkey.labels = c(
    expression("1"^" "),
    expression("10"^"-1"),
    expression("10"^"-2"),
    expression("10"^"-8")
  ), 
  description = 'Dotmap created by BoutrosLab.plotting.general',
  resolution = 50
)

dev.off()



# Writing a results table -------------------------------------------------

setdiff(peak.variants$ID, all.results$ID)
# [1] "rs2941509:chr17:39764941:T:C"

save(sig.results, all.results, file = "Database.m6A.analysis/peak.variants.allelic.imbalance.rsav")
