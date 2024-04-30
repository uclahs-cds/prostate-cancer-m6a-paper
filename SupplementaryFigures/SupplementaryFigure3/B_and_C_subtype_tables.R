library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

source('plot.general.contingency.table.R');

cna.clusters <- read.delim(
	file = '2022-01-18_CPCGENE_cluster_assignment.txt',
	as.is = TRUE
	);
cna <- cna.clusters[, c('sample', 'original')];

rna <- read.delim(
	file = '2018-07-13_rnaseq_consensus_cluster_assignment.txt',
	as.is = TRUE
	);

# load m6A clusters
m6a <- read.delim(
	file = '2022-01-18_m6a_adjusted_peak_cluster_assignment.txt',
	as.is = TRUE
	);

# overlap of 92 samples with RNA/protein
rna$ID <- paste0(rna$ID, '-F1');
patient.overlap.rna <- intersect(rna$ID, m6a$sample);
rna.clusters <- rna[match(patient.overlap.rna, rna$ID), ];
m6a.clusters.rna <- m6a[match(patient.overlap.rna, m6a$sample), ];

patient.overlap.cna <- intersect(cna$sample, m6a$sample);
cna.clusters <- cna[match(patient.overlap.cna, cna$sample), ];
m6a.clusters.cna <- m6a[match(patient.overlap.cna, m6a$sample), ];

plot.data.cna <- table(
	paste0('P', m6a.clusters.cna$cluster),
	paste0('C', cna.clusters$original)
	);

# remove C1, C2 because there is less than 10 samples that overlap (effects statistics)
cna.clusters.filtered <- cna.clusters[-which(cna.clusters$original == 1 | cna.clusters$original == 2), ];
m6a.clusters.filtered <- m6a.clusters.cna[match(cna.clusters.filtered$sample, m6a.clusters.cna$sample), ];
plot.data.cna <- table(
	paste0('P', m6a.clusters.filtered$cluster),
	paste0('C', cna.clusters.filtered$original)
	);

plot.general.contingency.table(
	plot.data = plot.data.cna,
	filename = 'SupplementaryFigure3C.png',
	xlabel = 'CNA subtype',
	ylabel = expression(bold('m')^bold('6') * bold('A subtype'))
	);

plot.data.rna <- table(
	paste0('P', m6a.clusters.rna$cluster),
	paste0('R', rna.clusters$RNA_cluster)
	);
plot.general.contingency.table(
	plot.data = plot.data.rna,
	filename = 'SupplementaryFigure3B.png',
	xlabel = 'RNA subtype',
	ylabel = expression(bold('m')^bold('6') * bold('A subtype'))
	);
