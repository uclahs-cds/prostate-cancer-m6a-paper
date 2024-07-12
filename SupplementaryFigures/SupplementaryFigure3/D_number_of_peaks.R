library(BoutrosLab.plotting.general);

# read in significant univariate results
s.peaks <- read.delim(
    file = '2022-04-29_m6A_peaks_with_significant_clinical_associations.txt',
    as.is = TRUE
    );

covs.to.keep <- c('pISUP', 'PSA', 'pT', 'age', 'IDC');

#count peaks associated with clinical covariates
cov.count <- table(s.peaks$clinical.feature);

for (cov in rev(covs.to.keep)) {
    indices <- grepl(cov, names(cov.count), ignore.case = TRUE);
    if (any(indices)) {
        names(cov.count)[which(indices)] <- cov;
        }
    }
cov.count <- cov.count[covs.to.keep];
cov.count[is.na(cov.count)] <- 0;

#match the names of the clinical covariates to those in the cluster-focused barplot
names(cov.count) <- covs.to.keep
names(cov.count) <- sub('age', 'Age', names(cov.count));
names(cov.count) <- sub('IDC', 'IDC/CA', names(cov.count));

cov.count.ordered <- cov.count[order(cov.count)]

toplot <- data.frame(
    count = as.numeric(cov.count.ordered),
    order = 1:length(covs.to.keep),
    stringsAsFactors = FALSE
    );

create.barplot(
    formula = order ~ count,
    data = toplot,
    filename = 'SupplementaryFigure3D.png',
    xaxis.tck = 0.5,
    yaxis.tck = 0.5,
    xlab.cex = 0.7,
    xlab.label = 'Peaks\n(Q < 0.1)',
    ylab.cex = 0,
    yaxis.lab = names(cov.count.ordered),
    yaxis.rot = 0,
    xaxis.cex = 0.6,
    yaxis.cex = 0.6,
    xlimit = c(0, 15),
    xat = seq(0, 15, 5),
    bottom.padding = 0,
    top.padding = 0,
    plot.horizontal = TRUE
    );
