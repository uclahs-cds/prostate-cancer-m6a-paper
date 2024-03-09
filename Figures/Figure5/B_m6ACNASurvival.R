# Analysis and visualization relating to Figure 5B

library(prostate.cancer.CNA.datasets);
library(BoutrosLab.plotting.general);
library(meta);
library(survival);

### path for saving figures and data needed to make table
path.out <- '/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/analysis/GRCh38/CNA_BCR/';

### path to m6A_genes.csv
path.m6a <- '/hot/project/disease/ProstateTumor/PRAD-000080-N6Methyl/m6A/analysis/GRCh38/CNA_BCR/m6A_genes.csv';

m6A <- read.csv(path.m6a, header = TRUE)
get.m6A.rate <- function(data.CNA) {
    # cohort name
    data.name <- deparse(substitute(data.CNA));
    name <- gsub('\\..*', '', data.name);
    
    # only calculate mutation rate among patients with non-missing BCR
    clinical <- get(paste0(name, '.clinical'));
    clinical <- clinical[!is.na(clinical$bcr.bin),];
    data.cna.clinical <- merge(data.CNA, clinical, by = 'patient.id');

    m6A.rate <- m6A;
    m6A.rate$rate <- NA;
    for(i in 1:nrow(m6A.rate)) {
        # some of the cohorts name these genes differently
        gene.name <- m6A.rate[i, 'symbol'];
        gene.name <- ifelse(
            test = gene.name %in% c('KIAA1429', 'VIRMA'),
            yes = 'KIAA1429|VIRMA',
            no = ifelse(
                test = gene.name %in% c('METTL16', 'METT10D'),
                yes = 'METTL16|METT10D',
                no = gene.name
                )
            );

        gene.index <- grep(gene.name, colnames(data.cna.clinical));
        # hieronymus has 2 duplicated HNRNPA2B1 cols
        gene.index <- ifelse(name == 'hieronymus' & gene.name == 'HNRNPA2B1', gene.index[1], gene.index);
        stopifnot(length(gene.index) == 1);
        
        if (m6A.rate[i, 'mutation.type'] == 'Gain') {
            m6A.rate$rate[i] <- mean(data.cna.clinical[, gene.index] == 1, na.rm = TRUE);
            }
        else {
            m6A.rate$rate[i] <- mean(data.cna.clinical[, gene.index] == -1, na.rm = TRUE);
            }
        }
    m6A.rate$cohort <- rep(name, nrow(m6A.rate));
    
    # Reported sample size should be number of nonmissing BCR and CNA
    m6A.rate$cohort[m6A.rate$cohort == 'cpcgene'] <- paste('ICGC PRAD-CA n=', sum(!is.na(analysis.cpcgene$analysis$bcr.bin)));
    m6A.rate$cohort[m6A.rate$cohort == 'taylor2010'] <- paste('Taylor et al. n=', sum(!is.na(analysis.taylor$analysis$bcr.bin)));
    m6A.rate$cohort[m6A.rate$cohort == 'hieronymus2014'] <- paste('Hieronymus et al. n=', sum(!is.na(analysis.hieronymus2014$analysis$bcr.bin)));
    m6A.rate$cohort[m6A.rate$cohort == 'cambridge'] <- paste('Ross-Adams et al. #1 n=', sum(!is.na(analysis.cambridge$analysis$bcr.bin)));
    m6A.rate$cohort[m6A.rate$cohort == 'stockholm'] <- paste('Ross-Adams et al. #2 n=', sum(!is.na(analysis.stockholm$analysis$bcr.bin)));
    m6A.rate$cohort[m6A.rate$cohort == 'TCGA'] <- paste('TCGA PRAD n=', sum(!is.na(analysis.TCGA$analysis$bcr.bin)));
    return(m6A.rate);
    }
# TCGA
clinical.t <- merge(TCGA.clinical, TCGA.CNA, by = 'patient.id');
gene.t <- clinical.t[, c('patient.id', m6A$symbol)];
gene.t[m6A$symbol] <- lapply(
    X = gene.t[m6A$symbol],
    FUN = function(x) {
        factor(x, levels = c(0, -1, 1), labels = c('Neutral', 'Loss', 'Gain'))
        }
    );

# Calculate frequency
freq.raw.TCGA <- lapply(gene.t[m6A$symbol], table);
freq.TCGA <- do.call(rbind, freq.raw.TCGA);

bcr <- clinical.t[, c('patient.id', 'bcr.bin', 'bcr.month')];
analysis.TCGA.raw <- merge(bcr, gene.t, by = 'patient.id');
analysis.TCGA <- list(analysis = analysis.TCGA.raw, freq = freq.TCGA);

# cpcgene
m6A$symbol[which(m6A$symbol == 'VIRMA')] <- 'KIAA1429';
clinical <- merge(cpcgene.clinical, cpcgene.CNA, by = 'patient.id');
gene.f <- clinical[, c('patient.id', m6A$symbol)];
gene.f[m6A$symbol] <- lapply(
    X = gene.f[m6A$symbol],
    FUN = function(x) {
        factor(x, levels = c(0, -1, 1), labels = c('Neutral', 'Loss', 'Gain'))
        }
    );
# Calculate frequency
freq.raw.cpcg <- lapply(gene.f[m6A$symbol], table);
freq.cpcg <- do.call(rbind, freq.raw.cpcg);
colnames(gene.f)[which(colnames(gene.f) == 'KIAA1429')] <- 'VIRMA';
rownames(freq.cpcg)[which(rownames(freq.cpcg) == 'KIAA1429')] <- 'VIRMA';

bcr <- clinical[, c('patient.id', 'bcr.bin', 'bcr.month')];
analysis.cpcgene.raw <- merge(bcr, gene.f, by = 'patient.id');
analysis.cpcgene <- list(analysis = analysis.cpcgene.raw, freq = freq.cpcg);

# Other cohorts
m6A$symbol[which(m6A$symbol == 'METTL16')] <- 'METT10D';
## Generate dataset for model building
data.analysis <- function(data.CNA, data.clinical) {
    clinical <- merge(data.clinical, data.CNA, by = 'patient.id');
    gene <- clinical[, c('patient.id', m6A$symbol)];
    gene[m6A$symbol] <- lapply(
        X = gene[m6A$symbol],
        FUN = function(x) {
            factor(x, levels = c(0, -1, 1), labels = c('Neutral', 'Loss', 'Gain'))
            }
        );
    # Calculate frequency
    freq.raw <- lapply(gene[m6A$symbol], table);
    freq <- do.call(rbind, freq.raw);
    colnames(gene)[which(colnames(gene) == 'KIAA1429')] <- 'VIRMA';
    colnames(gene)[which(colnames(gene) == 'METT10D')] <- 'METTL16';
    rownames(freq)[which(rownames(freq) == 'KIAA1429')] <- 'VIRMA';
    rownames(freq)[which(rownames(freq) == 'METT10D')] <- 'METTL16';
    bcr <- clinical[, c('patient.id', 'bcr.bin', 'bcr.month')];
    analysis <- merge(bcr, gene, by = 'patient.id');
    return.list <- list(analysis = analysis, freq = freq);
    return(return.list);
    }
analysis.stockholm <- data.analysis(stockholm.CNA, stockholm.clinical);
analysis.cambridge <- data.analysis(cambridge.CNA, cambridge.clinical);
analysis.taylor <- data.analysis(taylor2010.CNA, taylor2010.clinical);
analysis.hieronymus2014 <- data.analysis(hieronymus2014.CNA, hieronymus2014.clinical);

stockholm.m6A <- get.m6A.rate(stockholm.CNA);
cambridge.m6A <- get.m6A.rate(cambridge.CNA);
hieronymus.m6A <- get.m6A.rate(hieronymus2014.CNA);
taylor.m6A <- get.m6A.rate(taylor2010.CNA);
TCGA.m6A <- get.m6A.rate(TCGA.CNA);
cpcgene.m6A <- get.m6A.rate(cpcgene.CNA);

rate.m6A <- rbind(cpcgene.m6A, TCGA.m6A, stockholm.m6A, cambridge.m6A, hieronymus.m6A, taylor.m6A);

# data for table
analysis.cpcgene$analysis$cohort <- 'ICGC PRAD-CA';
analysis.TCGA$analysis$cohort <- 'TCGA PRAD';
analysis.taylor$analysis$cohort <- 'Taylor et al.';
analysis.hieronymus2014$analysis$cohort <- 'Hieronymus et al.';
analysis.cambridge$analysis$cohort <- 'Ross-Adams et al. #1';
analysis.stockholm$analysis$cohort <- 'Ross-Adams et al. #2';

analysis <- rbind(analysis.cpcgene$analysis, analysis.TCGA$analysis, analysis.taylor$analysis, analysis.hieronymus2014$analysis, analysis.cambridge$analysis, analysis.stockholm$analysis);
analysis$BCR <- factor(analysis$bcr.bin, levels = c(1, 0), labels = c('Yes', 'No'))

# only keep patients with non-missing bcr
analysis <- analysis[!is.na(analysis$bcr.bin),];

save(analysis, file = file.path(path.out, 'combined_data.RData'));
##########

# add overall mutation rates to rate.m6A
rate.m6A$symbol[which(rate.m6A$symbol == 'KIAA1429')] <- 'VIRMA';
rate.m6A$symbol[which(rate.m6A$symbol == 'METT10D')] <- 'METTL16';
overall.rates <- lapply(
    X = unique(rate.m6A$symbol),
    FUN = function(x) {
        temp <- rate.m6A[which(rate.m6A$symbol == x)[1],];
        temp$cohort <- 'Overall';
        temp$rate <- mean(analysis[,x] == temp$mutation.type, na.rm = TRUE);
        return(temp);
        }
    );
overall.rates <- do.call(rbind, overall.rates);
rate.m6A <- rbind(rate.m6A, overall.rates);
rate.m6A$cohort[rate.m6A$cohort == 'Overall'] <- paste('Overall n=', nrow(analysis));

# Calculate hazard ratio
get.hr.data <- function(data, gene) {
    if (data$freq[gene, ]['Gain'] < 5 & data$freq[gene, ]['Loss'] < 5) {
        index <- c('coef', 'exp.coef.', 'se.coef.',  'z',  'Pr...z..', 'lower..95', 'upper..95');
        coef <- data.frame();
        coef[1:2, 1:length(index)] <- NA;
        colnames(coef) <- index;
        rownames(coef) <- c('Loss', 'Gain');
        }
    else if (data$freq[gene, ]['Gain'] < 5 & data$freq[gene, ]['Loss'] >= 5) {
        data.new <- data$analysis[which(data$analysis[, gene] != 'Gain'), ];
        mod <- coxph(Surv(bcr.month, bcr.bin) ~ get(gene) , data = data.new, x = TRUE);
        sum.mod <- summary(mod);
        coef <- data.frame(sum.mod$coefficients);
        int <- data.frame(sum.mod$conf.int);
        coef <- cbind(coef, int[, c('lower..95', 'upper..95')]);
        rownames(coef) <- gsub('get\\(gene\\)', '', rownames(coef));
        coef['Gain', ] <- rep(NA, ncol(coef));
        }
    else if( data$freq[gene, ]['Loss'] < 5 & data$freq[gene, ]['Gain'] >= 5) {
            data.new <- data$analysis[which(data$analysis[, gene] != 'Loss'), ]
            mod <- coxph(Surv(bcr.month, bcr.bin) ~ get(gene), data = data.new, x = TRUE);
            sum.mod <- summary(mod);
            coef <- data.frame(sum.mod$coefficients);
            int <- data.frame(sum.mod$conf.int);
            coef <- cbind(coef, int[, c('lower..95', 'upper..95')]);
            rownames(coef) <- gsub('get\\(gene\\)', '', rownames(coef));
            coef['Loss', ] <- rep(NA, ncol(coef));
            }
    else {
            mod <- coxph(Surv(bcr.month, bcr.bin) ~ get(gene), data = data$analysis, x = TRUE);
            sum.mod <- summary(mod);
            coef <- data.frame(sum.mod$coefficients);
            int <- data.frame(sum.mod$conf.int);
            coef <- cbind(coef, int[, c('lower..95', 'upper..95')]);
            rownames(coef) <- gsub('get\\(gene\\)', '', rownames(coef));
            }
    coef <- data.frame(lapply(coef, as.numeric), row.names = rownames(coef));
    data.name <- deparse(substitute(data));
    name <- gsub('.*\\.', '', data.name);
    coef$cohort <- rep(name, nrow(coef));
    coef$genenames <- rep(gene, nrow(coef));
    coef$type <- c('Loss', 'Gain');
    
    # fix bug due to genes with multiple names
     gene.name <- ifelse(
            test = gene %in% c('KIAA1429', 'VIRMA'),
            yes = 'KIAA1429|VIRMA',
            no = ifelse(
                test = gene %in% c('METTL16', 'METT10D'),
                yes = 'METTL16|METT10D',
                no = gene
                )
            );
    gene.index <- grep(gene.name, m6A$symbol);
    stopifnot(length(gene.index) == 1);
    coef$mutationtype <- rep(m6A[gene.index, 'mutation.type'], 2)
    
    stopifnot(rownames(coef) == c('Loss', 'Gain'));
    return(coef);
    }

#To test get.hr.data function
# get.hr.data(data = analysis.hieronymus2014, gene = 'VIRMA');

# data preparation
hr.gene <- function(gene) {
    hr.TCGA <- get.hr.data(analysis.TCGA, gene);
    hr.cpcg <- get.hr.data(analysis.cpcgene, gene);
    hr.stockholm <- get.hr.data(analysis.stockholm, gene);
    hr.cambridge <- get.hr.data(analysis.cambridge, gene);
    hr.hieronymus2014 <- get.hr.data(analysis.hieronymus2014, gene);
    hr.taylor <- get.hr.data(analysis.taylor, gene);
    hr <- rbind(hr.TCGA, hr.cpcg, hr.stockholm, hr.cambridge, hr.hieronymus2014, hr.taylor);
    ## Overall hazard ratio for HNRNPA2B1
    #hr.loss <- hr[seq(1, nrow(hr) - 1, 2), ];
    #hr.gain <- hr[seq(2, nrow(hr), 2), ];
    hr.loss <- hr[grepl('Loss', rownames(hr)), ];
    hr.gain <- hr[grepl('Gain', rownames(hr)), ];
    if (all(is.na(hr.loss$coef)) & any(!is.na(hr.gain$coef))) {
        hr.overall.gain <- metagen(TE = coef, seTE = se.coef., data = hr, subset = type == 'Gain', comb.fixed = FALSE);
        hr['OverallGain', ] <- c(hr.overall.gain$TE.random, exp(hr.overall.gain$TE.random), hr.overall.gain$seTE.random, hr.overall.gain$zval.random, hr.overall.gain$pval.random, exp(hr.overall.gain$lower.random), exp(hr.overall.gain$upper.random), 'Overall', gene, 'Gain', unique(hr$mutationtype));
        hr['OverallLoss', ] <- rep(NA, ncol(hr));
        hr['OverallLoss', c('cohort', 'genenames', 'type', 'mutationtype')] <- c('Overall', gene, 'Loss', na.omit(unique(hr$mutationtype)));
        }
    else if(any(!is.na(hr.loss$coef)) & all(is.na(hr.gain$coef))) {
        hr.overall.loss <- metagen(TE = coef, seTE = se.coef., data = hr, subset = type == 'Loss', comb.fixed = FALSE);
        hr['OverallLoss', ] <- c(hr.overall.loss$TE.random, exp(hr.overall.loss$TE.random), hr.overall.loss$seTE.random, hr.overall.loss$zval.random, hr.overall.loss$pval.random, exp(hr.overall.loss$lower.random), exp(hr.overall.loss$upper.random), 'Overall', gene, 'Loss', unique(hr$mutationtype));
        hr['OverallGain', ] <- rep(NA, ncol(hr));
        hr['OverallGain', c('cohort', 'genenames', 'type', 'mutationtype')] <- c('Overall', gene, 'Gain', na.omit(unique(hr$mutationtype)));
        }
    else{
        hr.overall.gain <- metagen(TE = coef, seTE = se.coef., data = hr, subset = type == 'Gain', comb.fixed = FALSE);
        hr.overall.loss <- metagen(TE = coef, seTE = se.coef., data = hr, subset = type == 'Loss', comb.fixed = FALSE);
        hr['OverallGain', ] <- c(hr.overall.gain$TE.random, exp(hr.overall.gain$TE.random), hr.overall.gain$seTE.random, hr.overall.gain$zval.random, hr.overall.gain$pval.random, exp(hr.overall.gain$lower.random), exp(hr.overall.gain$upper.random), 'Overall', gene, 'Gain', unique(hr$mutationtype));
        hr['OverallLoss', ] <- c(hr.overall.loss$TE.random, exp(hr.overall.loss$TE.random), hr.overall.loss$seTE.random, hr.overall.loss$zval.random, hr.overall.loss$pval.random, exp(hr.overall.loss$lower.random), exp(hr.overall.loss$upper.random), 'Overall', gene, 'Loss', unique(hr$mutationtype));
        }
    hr$coef <- as.numeric(hr$coef);
    hr$exp.coef. <- as.numeric(hr$exp.coef.);
    hr$lower..95 <- as.numeric(hr$lower..95);
    hr$upper..95 <- as.numeric(hr$upper..95);
    hr$coef.log2 <- log2(hr$exp.coef.);
    hr$upper <- log2(hr$upper..95);
    hr$lower <- log2(hr$lower..95);
    hr$logp <- -log10(as.numeric(hr$Pr...z..));
    
    # Reported sample size should be number of nonmissing BCR and CNA
    hr$cohort[hr$cohort == 'cpcgene'] <- paste('ICGC PRAD-CA n=', sum(!is.na(analysis.cpcgene$analysis$bcr.bin)));
    hr$cohort[hr$cohort == 'taylor'] <- paste('Taylor et al. n=', sum(!is.na(analysis.taylor$analysis$bcr.bin)));
    hr$cohort[hr$cohort == 'hieronymus2014'] <- paste('Hieronymus et al. n=', sum(!is.na(analysis.hieronymus2014$analysis$bcr.bin)));
    hr$cohort[hr$cohort == 'cambridge'] <- paste('Ross-Adams et al. #1 n=',  sum(!is.na(analysis.cambridge$analysis$bcr.bin)));
    hr$cohort[hr$cohort == 'stockholm'] <- paste('Ross-Adams et al. #2 n=', sum(!is.na(analysis.stockholm$analysis$bcr.bin)));
    hr$cohort[hr$cohort == 'TCGA'] <- paste('TCGA PRAD n=', sum(!is.na(analysis.TCGA$analysis$bcr.bin)));
    hr$cohort[hr$cohort == 'Overall'] <- paste('Overall n=', nrow(analysis));

    return(hr)
    }

cohort.labels <- c(
    as.expression(bquote(bold(.(paste0('TCGA PRAD n=', sum(!is.na(analysis.TCGA$analysis$bcr.bin))))))),
    as.expression(bquote(bold(.(paste0('ICGC PRAD-CA n=', sum(!is.na(analysis.cpcgene$analysis$bcr.bin))))))),
    as.expression(bquote(bold('Ross-Adams ' * bolditalic(et ~ al.) * ~ .(paste0('#2 n=', sum(!is.na(analysis.stockholm$analysis$bcr.bin)), ''))))),
    as.expression(bquote(bold('Ross-Adams ' * bolditalic(et ~ al.) * ~ .(paste0('#1 n=', sum(!is.na(analysis.cambridge$analysis$bcr.bin)), ''))))),
    as.expression(bquote(bold('Hieronymus ' * bolditalic(et ~ al.) * ~ .(paste0('n=', sum(!is.na(analysis.hieronymus2014$analysis$bcr.bin)), ''))))),
    as.expression(bquote(bold('Taylor ' * bolditalic(et ~ al.) * ~ .(paste0('n=', sum(!is.na(analysis.taylor$analysis$bcr.bin)), ''))))),
    as.expression(bquote(bold(.(paste0('Overall n=', nrow(analysis))))))
    );

# Plot Function
### testing: 
# ZC3H13.gene <- hr.gene('ZC3H13')
# create.plot.seg(gene = 'ZC3H13', hr = ZC3H13.gene);
# create.plot.seg(gene = 'HNRNPA2B1', hr = HNRNPA2B1.gene);
# create.plot.seg(gene = 'FTO', hr = FTO.gene);

create.plot.seg <- function(gene, hr) {
    hr$cohort <- factor(hr$cohort, levels = unique(hr$cohort));
    hr.loss <- hr[hr$type == 'Loss', ];
    hr.gain <- hr[hr$type == 'Gain', ];
    if(na.omit(unique(hr$mutationtype) == 'Gain')) {
        ## barplot for p-value
        if(ceiling(max(hr.gain$logp, na.rm = TRUE)) < 2) {
            pvalue.at <- seq(0, ceiling(max(hr.gain$logp, na.rm = TRUE)), by = 1);
            #if (all(pvalue.at) %in% c(0,1)) {
            #    pvalue.at <- c(0, 1, 2);
            #    }
            pvalue.at.labels <- sapply(
                X = pvalue.at,
                FUN = function(x) {
                    as.expression(bquote(bold('10'^-.(as.character(x)))))
                    }
                );
            }
        else {
            pvalue.at <- seq(0, ceiling(max(hr.gain$logp, na.rm = TRUE)), by = 2);
            #if (all(pvalue.at) %in% c(0,1)) {
            #    pvalue.at <- c(0, 1, 2);
            #    }
            pvalue.at.labels <- sapply(
                X = pvalue.at,
                FUN = function(x) {
                    as.expression(bquote(bold('10'^-.(as.character(x)))))
                    }
                );
            }
        barplot <- BoutrosLab.plotting.general::create.barplot(
            formula = as.formula(cohort ~ logp),
            data = hr.gain,
            disable.factor.sorting = TRUE,
            plot.horizontal = TRUE,
            xlab.label = 'P',
            ylab.label = '',
            xlimits = c(0, ceiling(max(hr.gain$logp, na.rm = TRUE))),
            xat = pvalue.at,
            xaxis.lab = pvalue.at.labels,
            xaxis.cex = .3,
            xlab.cex = 0.35,
            yaxis.lab = rep('', nrow(hr.gain)),
            yaxis.cex = 0,
            abline.v = -log10(0.05),
            abline.lty = 3,
            abline.col = 'grey',
            xaxis.tck = c(.5, 0),
            yaxis.tck = 0
            );
        data <- rate.m6A[rate.m6A$symbol == gene, ];
        data$cohort <- factor(data$cohort, levels = levels(hr$cohort));
        data <- data[order(data$cohort),];
        stopifnot(all(data$cohort == hr.gain$cohort));
        barplot.rate <- BoutrosLab.plotting.general::create.barplot(
            formula = as.formula(cohort ~ rate),
            data = data,
            disable.factor.sorting = TRUE,
            plot.horizontal = TRUE,
            xlab.label = 'Mutation Rate',
            xlimits = c(0, max(na.omit(data$rate)) + 0.1),
            xat = seq(0, 1, by = 0.1),
            xaxis.cex = 0.3,
            yaxis.cex = 0,
            xlab.cex = 0.35,
            ylab.cex = 0,
            abline.lty = 3,
            abline.col = 'grey',
            xaxis.tck = c(.5, 0),
            yaxis.tck = 0
            );
        ## Create segplot
        segplot <- BoutrosLab.plotting.general::create.segplot(
            formula = as.formula(cohort ~ lower + upper),
            data = hr.gain,
            resolution = 1000,
            abline.v = 0,
            abline.lty = 3,
            abline.col = 'grey',
            centers = hr.gain$coef.log2,
            xaxis.rot = 0,
            xaxis.cex = 0.4,
            yaxis.cex = 0.33,
            xlab.cex = 0.5,
            ylab.cex = 0,
            xlimits = c(-2, 4.1),
            xat = seq(-2, 4, 1),
            xaxis.lab = sapply(
                X = -2:4,
                FUN = function(x) {
                    as.expression(bquote(bold('2'^.(as.character(x)))))
                    }
                ),
            xlab.label = 'Hazard Ratio',
            yaxis.lab = cohort.labels,
            ylab.label = '',
            xaxis.tck = c(.5, 0),
            yaxis.tck = 0,
            symbol.cex = .65
            );
        create.multipanelplot(
            filename = file.path(path.out, paste0(Sys.Date(),'_', gene, '_segplot.tiff')),
            plot.objects = list(segplot, barplot, barplot.rate),
            layout.width = 3,
            layout.height = 1,
            plot.objects.widths = c(4, 1.1, 1.1),
            x.spacing = -3,
            left.padding = -0.7,
            width = 3.5,
            height = 2.5,
            left.legend.padding = 0,
            main = paste(gene, ' Gain'),
            main.cex = .6,
            top.padding = -2,
            resolution = 1000,
            main.y = -0.3,
            right.padding = -2.5,
            bottom.padding = -2
            );
        }
    else {
        if(ceiling(max(hr.loss$logp, na.rm = TRUE)) < 2) {
            pvalue.at <- seq(0, ceiling(max(hr.loss$logp, na.rm = TRUE)), by = 1);
            #if (all(pvalue.at) %in% c(0,1)) {
            #    pvalue.at <- c(0, 1, 2);
            #    }
            pvalue.at.labels <- sapply(
                X = pvalue.at,
                FUN = function(x) {
                    as.expression(bquote(bold('10'^-.(as.character(x)))))
                    }
                );
            }
        else {
            pvalue.at <- seq(0, ceiling(max(hr.loss$logp, na.rm = TRUE)), by = 2);
            #if (all(pvalue.at) %in% c(0,1)) {
            #    pvalue.at <- c(0, 1, 2);
            #    }
            pvalue.at.labels <- sapply(
                X = pvalue.at,
                FUN = function(x) {
                    as.expression(bquote(bold('10'^-.(as.character(x)))))
                    }
                );
            }
        barplot <- BoutrosLab.plotting.general::create.barplot(
            formula = as.formula(cohort ~ logp),
            data = hr.loss,
            disable.factor.sorting = TRUE,
            plot.horizontal = TRUE,
            xlab.label = 'P',
            ylab.label = '',
            xlimits = c(0, ceiling(max(hr.loss$logp, na.rm =TRUE))),
            xat = pvalue.at,
            xaxis.lab = pvalue.at.labels,
            xaxis.cex = .3,
            xlab.cex = 0.35,
            yaxis.lab = rep('', nrow(hr.loss)),
            yaxis.cex = 0,
            abline.v = -log10(0.05),
            abline.lty = 3,
            abline.col = 'grey',
            xaxis.tck = c(.5, 0),
            yaxis.tck = 0
            );
        data <- rate.m6A[rate.m6A$symbol == gene, ];
        data$cohort <- factor(data$cohort, levels = levels(hr$cohort));
        data <- data[order(data$cohort),];
        stopifnot(all(data$cohort == hr.loss$cohort));
        
        barplot.rate <- BoutrosLab.plotting.general::create.barplot(
            formula = as.formula(cohort ~ rate),
            data = data,
            disable.factor.sorting = TRUE,
            plot.horizontal = TRUE,
            xlab.label = 'Mutation Rate',
            xlimits = c(0, max(data$rate) + 0.1),
            xat = seq(0, 1, by = 0.1),
            xaxis.cex = 0.3,
            xlab.cex = 0.35,
            abline.lty = 3,
            abline.col = 'grey',
            ylab.cex = 0,
            yaxis.cex = 0,
            xaxis.tck = c(.5, 0),
            yaxis.tck = 0
            );
        segplot <- BoutrosLab.plotting.general::create.segplot(
            formula = as.formula(cohort ~ lower + upper),
            data = hr.loss,
            resolution = 1000,
            abline.v = 0,
            abline.lty = 3,
            abline.col = 'grey',
            centers = hr.loss$coef.log2,
            xaxis.rot = 0,
            xaxis.cex = 0.4,
            yaxis.cex = 0.33,
            xlab.cex = 0.5,
            ylab.cex = 0,
            xlimits = c(-2, 4.1),
            xat = seq(-2, 4, 1),
            xaxis.lab = sapply(
                X = -2:4,
                FUN = function(x) {
                    as.expression(bquote(bold('2'^.(as.character(x)))))
                    }
                ),
            xlab.label = 'Hazard Ratio',
            yaxis.lab = cohort.labels,
            ylab.label = '',
            xaxis.tck = c(.5, 0),
            yaxis.tck = 0,
            symbol.cex = .65
            );
        # Combine 2 plots
        create.multipanelplot(
            filename = file.path(path.out, paste0(Sys.Date(), '_', gene, '_segplot.tiff')),
            plot.objects = list(segplot, barplot, barplot.rate),
            layout.width = 3,
            layout.height = 1,
            plot.objects.widths = c(4, 1.1, 1.1),
            x.spacing = -3,
            left.padding = -0.7,
            width = 3.5,
            height = 2.5,
            left.legend.padding = 0,
            main = paste(gene, ' Loss'),
            main.cex = .6,
            top.padding = -2,
            resolution = 1000,
            main.y = -0.3,
            right.padding = -2.5,
            bottom.padding = -2
            );
        }
    }

HNRNPA2B1.gene <- hr.gene('HNRNPA2B1');
create.plot.seg('HNRNPA2B1', HNRNPA2B1.gene);

IGF2BP3.gene <- hr.gene('IGF2BP3');
IGF2BP3.segplot <- create.plot.seg('IGF2BP3', IGF2BP3.gene);

YTHDF3.gene <- hr.gene('YTHDF3');
create.plot.seg('YTHDF3', YTHDF3.gene);

VIRMA.gene <- hr.gene('VIRMA');
create.plot.seg('VIRMA', VIRMA.gene);

YTHDF1.gene <- hr.gene('YTHDF1');
create.plot.seg('YTHDF1', YTHDF1.gene);

YTHDC2.gene <- hr.gene('YTHDC2')
create.plot.seg('YTHDC2', YTHDC2.gene);

METTL16.gene <- hr.gene('METTL16')
create.plot.seg('METTL16', METTL16.gene);

FTO.gene <- hr.gene('FTO')
create.plot.seg('FTO', FTO.gene);

ZC3H13.gene <- hr.gene('ZC3H13')
create.plot.seg('ZC3H13', ZC3H13.gene);


################# Make forest plot for overall meta-analysis HRs
subset.data <- function(data) {
    if(na.omit(unique(data$mutationtype)) == 'Gain') {
        sub.data <- data['OverallGain', ]
        }
    else {
        sub.data <- data['OverallLoss', ]
        }
    }
HNRNPA2B1.gene.overall <- subset.data(HNRNPA2B1.gene);
IGF2BP3.gene.overall <- subset.data(IGF2BP3.gene);
YTHDF3.gene.overall <- subset.data(YTHDF3.gene);
VIRMA.gene.overall <- subset.data(VIRMA.gene);
YTHDF1.gene.overall <- subset.data(YTHDF1.gene);
YTHDC2.gene.overall <- subset.data(YTHDC2.gene);
METTL16.gene.overall <- subset.data(METTL16.gene);
FTO.gene.overall <- subset.data(FTO.gene)
ZC3H13.gene.overall <- subset.data(ZC3H13.gene)


##### forestplot for overall HR
overall.data <- rbind(
    HNRNPA2B1.gene.overall,
    IGF2BP3.gene.overall,
    YTHDF3.gene.overall,
    VIRMA.gene.overall,
    YTHDF1.gene.overall,
    YTHDC2.gene.overall,
    METTL16.gene.overall,
    FTO.gene.overall,
    ZC3H13.gene.overall
    );
overall.data <- overall.data[order(overall.data$coef), ];
overall.data$genenames <- factor(overall.data$genenames, levels = overall.data$genenames);
save(overall.data, file = file.path(path.out, 'overall_data.RData'));

pvalue.at <- seq(0, ceiling(max(overall.data$logp)), by = 2);
pvalue.at.labels <- sapply(
    X = pvalue.at,
    FUN = function(x) {
        as.expression(bquote(bold('10'^-.(as.character(x)))))
        }
    );
barplot <- BoutrosLab.plotting.general::create.barplot(
    formula = as.formula(genenames ~ logp),
    data = overall.data,
    disable.factor.sorting = TRUE,
    plot.horizontal = TRUE,
    xlab.label = 'P',
    ylab.label = '',
    xlimits = c(0, max(overall.data$logp) + 0.1),
    xat = pvalue.at,
    xaxis.lab = pvalue.at.labels,
    xaxis.cex = 0.365,
    xlab.cex = 0.45,
    yaxis.lab = rep('', nrow(overall.data)),
    yaxis.cex = 0,
    abline.v = -log10(0.05),
    abline.lty = 3,
    abline.col = 'grey',
    xaxis.tck = c(.5, 0),
    yaxis.tck = 0
    );
segplot <- BoutrosLab.plotting.general::create.segplot(
    formula = as.formula(genenames ~ lower + upper),
    data = overall.data,
    resolution = 50,
    abline.v = 0,
    abline.lty = 3,
    abline.col = 'grey',
    centers = overall.data$coef.log2,
    xaxis.rot = 0,
    xlab.cex = .45,
    xaxis.cex = .45,
    yaxis.cex = .4,
    yaxis.tck = 0,
    xaxis.tck = c(.5, 0),
    ylab.cex = 0,
    xlimits = c(-1, 2.1),
    xat = seq(-1, 2, 1),
    xaxis.lab = sapply(
        X = -1:2,
        FUN = function(x) {
            as.expression(bquote(bold('2'^.(as.character(x)))))
            }
        ),
    yaxis.lab = rep('', nrow(overall.data)),
    xlab.label = 'Hazard Ratio',
    ylab.label = 'Gene',
    symbol.cex = 0.7
    );
overall.rates <- rate.m6A[grepl('^Overall', rate.m6A$cohort),];
overall.rates$symbol <- factor(overall.rates$symbol, levels = levels(overall.data$genenames));
overall.rates <- overall.rates[order(overall.rates$symbol),];

stopifnot(all(overall.rates$symbol == overall.data$genenames));

overall.rates$ensyme[overall.rates$symbol %in% c('HNRNPA2B1', 'IGF2BP3', 'YTHDF3', 'YTHDF1', 'YTHDC2')] <- 'Reader';
overall.rates$ensyme[overall.rates$symbol %in% c('VIRMA', 'METTL16', 'ZC3H13')] <- 'Writer';
overall.rates$ensyme[overall.rates$symbol == 'FTO'] <- 'Eraser';

col <- c('#ff3030','#1e90ff')[match(overall.rates$mutation.type, c('Gain', 'Loss'))];
feature.col <- c('#ffa500', '#458b00', '#68228b')[match(overall.rates$ensyme, c('Eraser', 'Reader', 'Writer'))];
col <- cbind(col, feature.col);
col <- data.frame(col);
rownames(col) <- overall.rates$symbol;
total.colours <- 6;


covariates <- create.heatmap(
    t(col),
    clustering.method = 'none',
    total.colours = total.colours,
    print.colour.key = FALSE,
    yaxis.tck = 0,
    xaxis.tck = 0,
    input.colours = TRUE,
    yaxis.lab = NA,
    yaxis.cex = 0.4,
    axes.lwd = 0.8
    );
covariates

barplot.rate <- BoutrosLab.plotting.general::create.barplot(
    formula = as.formula(symbol ~ rate),
    data = overall.rates,
    disable.factor.sorting = TRUE,
    plot.horizontal = TRUE,
    xlab.label = 'Mutation Rate',
    xlimits = c(0, max(na.omit(overall.rates$rate)) + 0.1),
    xat = seq(0, 1, by = 0.1),
    xaxis.cex = 0.365,
    yaxis.cex = 0,
    xlab.cex = 0.45,
    ylab.cex = 0,
    abline.lty = 3,
    abline.col = 'grey',
    xaxis.tck = c(.5, 0),
    yaxis.tck = 0
    );

cov.legend <- list(
    legend = list(
        colours = c('#ff3030',"#1e90ff"),
        labels = c('Gain', 'Loss'),
        title = expression(underline(bold('Mutations'))),
        cex = .3
        ),
    legend = list(
        colours = c("#ffa500", "#458b00", "#68228b"),
        labels = c('Eraser', 'Reader', "Writer"),
        title = expression(underline(bold(m^'6'*A*paste(' ', 'Enzyme')))),
        cex = .3
        )
    );
legend.grob <- legend.grob(
    legends = cov.legend,
    title.just = 'left',
    label.cex = 0.3,
    # layout = c(2,1),
    title.cex = 0.3,
    size = .9
    );
plot.objects.widths <- c(1, 2, 1.2, 1.2);
# Combine 3 plots
create.multipanelplot(
    filename = file.path(path.out, paste0(Sys.Date(),'_meta_analysis_segplot.tiff')),
    plot.objects = list(covariates, segplot, barplot, barplot.rate),
    layout.width = 4,
    layout.height = 1,
    # layout.skip =  c(FALSE, FALSE, FALSE, FALSE),
    plot.objects.widths = plot.objects.widths,
    resolution = 300,
    legend = list(
        left = list(
            fun = legend.grob
            )
        ),
    left.padding = 0,
    right.legend.padding = 0,
    top.legend.padding = 0,
    height = 3.1,
    width = 4.2,
    x.spacing = -2,
    left.legend.padding = 0,
    top.padding = -1,
    bottom.padding = -1.5
    );