library(GenomicRanges);
library(BoutrosLab.utilities);
library(valr);

config <- 'D_E.config';
config <- read.config.file(config);

gene_id <- config$gene
gtf <- read_gtf(config$gtf, zero_based = FALSE);

# Creating a GRanges Object
genes.gr <- gtf[grepl(gene_id, gtf$gene_id) & 
                gtf$type == 'gene',
                c('chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name')];
genes.gr <- makeGRangesFromDataFrame(genes.gr, keep.extra.columns = T);

rm(gtf);

sample.paths <- unique(sub(
    '_I.+.bigWig',
    '',
    Sys.glob(paste0(config$patient.bigWigs, '*.bigWig'))
    ));
sample.paths <- sample.paths[grep('CPCG0', sample.paths)];
names(sample.paths) <- basename(sample.paths);

cell.line.paths <- unique(sub(
    '_I.+.bigWig',
    '',
    grep('PC-3|V16A', grep('-normoxia|-hypoxia', Sys.glob(paste0(config$cell.line.bigWigs, '*.bigWig')), value = T), value = T)
    ));
names(cell.line.paths) <- basename(cell.line.paths);
# Calculate coverage
calculate.coverage <- function(f.cov, r.cov, genes.gr) {
    histogram.coverage <- vector('list', length(genes.gr))
    names(histogram.coverage) <- genes.gr$gene_id
    for (x in seq_along(genes.gr)) {
        gr = genes.gr[x]
        bins = unlist(GenomicRanges::tile(x = gr, width = 1))
        GenomeInfoDb::seqlevels(bins) <- GenomeInfoDb::seqlevels(f.cov)
        if (as.character(strand(gr)) == '+') {
            cvg <- binnedAverage(
                bins = bins,
                numvar = f.cov,
                varname = 'cvg'
                )
            } else {
            cvg <- binnedAverage(
                bins = bins,
                numvar = r.cov,
                varname = 'cvg'
                )
            }
        histogram.coverage[[gr$gene_id]] <- cvg$cvg
        }
    return(histogram.coverage)
    }

# https://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r
slidingWindowAverage <- function(data, window, step){
    total <- length(data)
    spots <- seq(from=1, to=(total-window), by=step)
    result <- vector(length = length(spots))
    for(i in 1:length(spots)){
        result[i] <- mean(data[spots[i]:(spots[i]+window)])
        }
    return(result)
    }

process.bigWigs <- function(sample.path, type, genes.gr, flip = FALSE) {
    cat(sample.path, '\n')
    # Forward BigWig
    f.bw.file <- paste0(sample.path, '_', type, '.forward.bigWig');
    f.bw <- valr::read_bigwig(f.bw.file, set_strand = '+');
    f.gr <- makeGRangesFromDataFrame(f.bw, keep.extra.columns = T);
    f.gr <- f.gr[f.gr$score > 0]
    f.cov <- coverage(f.gr, weight = 'score')

    # Reverse BigWig
    r.bw.file <- paste0(sample.path, '_', type, '.reverse.bigWig')
    r.bw <- valr::read_bigwig(r.bw.file, set_strand = '-')
    r.gr <- makeGRangesFromDataFrame(r.bw, keep.extra.columns = T)
    r.gr <- r.gr[r.gr$score > 0]
    r.cov <- coverage(r.gr, weight = 'score')

    # Calculate coverage
    if (flip) { #strand for cell line BigWigs needs to be flipped due to how it was generated with deepTools
        coverage <- calculate.coverage(r.cov, f.cov, genes.gr)
        } else {
        coverage <- calculate.coverage(f.cov, r.cov, genes.gr)
        }
    smooth.data <- sapply(coverage, function(x) slidingWindowAverage(c(0, 0, 0, x, 0, 0, 0), 6, 1));
    return(smooth.data);
    }

# Input Bigwigs
input.coverage <- sapply(sample.paths, process.bigWigs, genes.gr = genes.gr, type = 'Input');
IP.coverage <- sapply(sample.paths, process.bigWigs, genes.gr = genes.gr, type = 'IP');
cell.input.coverage <- sapply(cell.line.paths, process.bigWigs, genes.gr = genes.gr, type = 'Input', flip = TRUE);
cell.IP.coverage <- sapply(cell.line.paths, process.bigWigs, genes.gr = genes.gr, type = 'IP', flip = TRUE);
colnames(input.coverage) <- names(sample.paths);
colnames(IP.coverage) <- names(sample.paths);
colnames(cell.input.coverage) <- names(cell.line.paths);
colnames(cell.IP.coverage) <- names(cell.line.paths);

save(IP.coverage, input.coverage, cell.input.coverage, cell.IP.coverage,
    file = paste0(gene_id,'_smooth_coverage.rsav')
    );

