sample.data <- read.delim('SACseq_m6Aseq_sample_overlap_unadjusted.txt')[-2, ]

create.scatterplot(
    SAC.seq ~ m6A.seq, 
    data = sample.data,
    filename = 'm6A_peak_vs_SAC_site_number_cor.pdf',
    cex = .75,
    xlab.label = expression(bold("m"^"6"*"A-seq"~"peaks")),
    ylab.label = 'SAC-seq sites',
    ylab.cex = .7,
    yaxis.cex = .45,
    xlab.cex = .7,
    xaxis.cex = .45,
    xlimits = c(-500, 12000),
    xat = seq(0, 12000, 2000),
    yat = seq(0, 12000, 2000),
    ylimits = c(-500, 12000),
    main.cex = .7,
    alpha = .55,
    xaxis.tck = .2,
    yaxis.tck = .2,
    type = c('p'),
    #add correlation
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
            key = get.corr.key(
                # two vectors of the same length
                y = sample.data$SAC.seq,
                x = sample.data$m6A.seq,
                # specify what to include in the key
                label.items = c('pearson','pearson.p'),
                # format key
                alpha.background = 0,
                key.cex = .6,
                key.title = "Pearson's correlation",
                title.cex = .6
                )
            ),
            # place key
            key.x = 0.1,
            key.y = 1,
            corner = c(0,1)
            )
        ),
    resolution = 300,
    width = 2.75,
    height = 2.5
    );
