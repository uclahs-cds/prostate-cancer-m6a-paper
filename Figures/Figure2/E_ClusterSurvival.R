#Perform analysis for and generate Figure 2E 
library(BoutrosLab.plotting.general);
library(BoutrosLab.plotting.survival);
library(BoutrosLab.utilities);

level <- 'peak'
input <- 'adjusted'
sample.k <- 5

#load m6A sample clusters (Figure 2A)
sample.part <- paste0(input, '_', level, '_sample_ConsensusClusterPlus');
sample.file <- paste0('/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/clustering/', sample.part, '/', sample.part, '.k=', sample.k, '.consensusClass.csv');
sample.clusters <- read.delim(
    file = sample.file,
    as.is = TRUE,
    header = FALSE,
    sep = ','
    );

#patient clinical data
clinical <- read.delim(
    file = '/hot/ref/database/ProstateTumor/Boutros-Yamaguchi-PRAD-CPCG/clinical/2020-02-25_CPC-GENE_clinical_from_database.txt',
    as.is = TRUE
    );

annot <- clinical[match(sample.clusters$V1, clinical$patient_id), ];
annot$cluster <- paste0('P', sample.clusters$V2);

#perform survival analysis on patient m6A clusters
annot.survival <- annot[-which(is.na(annot$time_to_bcr) | is.na(annot$bcr)), c('patient_id', 'bcr', 'time_to_bcr', 'cluster')];
surv <- Surv(
    time = annot.survival$time_to_bcr,
    event = as.numeric(annot.survival$bcr)
    );

#Run Cox proportional hazards model
fit <- coxph(
    formula = Surv(time_to_bcr, bcr) ~ cluster,
    data = annot.survival
    );
fit.summary <- summary(fit);

#Create Kaplan-Meier plot
cluster.km <- create.km.plot(
    survival.object = surv,
    patient.groups = annot.survival$cluster,
    xlab.label = 'Time (Years)',
    key.groups.title.cex = .8,
    key.groups.cex = .8,
    key.stats.cex = .8,
    left.padding = 5,
    censoring.pch.cex = .8,
    lwd = .9,
    ylab.label = 'Biochemical Relapse-Free Rate',
    xlab.cex = .8,
    ylab.cex = .8,
    xaxis.cex = .8,
    yaxis.cex = .8,
    show.risktable = TRUE,
    show.key.groups = TRUE,
    risktable.fontsize = 9,
    statistical.method = 'none'
    );

cluster.km$x.scales$tck <- c(.5, 0);
cluster.km$y.scales$tck <- c(.5, 0);

# Start custom output for km plot
png(
    file = generate.filename('m6A', 'cluster_survival_km_plot', 'png'),
    height = 4.5,
    width = 4.5,
    units = 'in',
    res = 300
    );

par(
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0),
    bg = 'white',
    cex = 1.5
    );

#Add custom text
plot.new();
plot(cluster.km, newpage = FALSE);

group.names <- c('clusterP2', 'clusterP3', 'clusterP4', 'clusterP5');
y <- 0.525;
for (i in 1:length(group.names)) {
    text(
        x = 0.475,
        y = y,
        adj = c(0, 1),
        cex = 0.5,
        labels = paste0(
            'HR: ', sprintf("%.2f", round(fit.summary$conf.int[group.names[i], 'exp(coef)'], digits = 2)),
            ' (', sprintf("%.2f", round(fit.summary$conf.int[group.names[i], 'lower .95'], digits = 2)),
            ', ', sprintf("%.2f", round(fit.summary$conf.int[group.names[i], 'upper .95'], digits = 2)), ') ',
            'P: ', sprintf("%.2f", round(fit.summary$coefficients[, 'Pr(>|z|)'][i], digits = 2))
            )
        );
    y <- y - 0.03;
    }

dev.off();






