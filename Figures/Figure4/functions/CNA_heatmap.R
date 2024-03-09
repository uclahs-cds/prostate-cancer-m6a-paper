#convert data in the format gene.loss gene.gain in separate columns into 1, -1, with 0 == neutral
# TO DO in future
make.tertiary.CNA <- function() {
    mut.mat.template <- matrix(0, nrow = nrow(mutation.data), ncol = length(unique(mut.labels)));
    colnames(mut.mat.template) <- unique(mut.labels);
    row.names(mut.mat.template) <- row.names(mutation.data);

    mut.list <- lapply(
        unique(mutation.types), function(mutation) {
            mut.mat <- mut.mat.template;
            mut.mat[, mut.labels[mutation.types == mutation]] <- mutation.data[, mutation.types == mutation];
            return(t(mut.mat));
            }
        );
    names(mut.list) <- unique(mutation.types);
    mut.list$CNA <- ifelse(mut.list$Loss == 1, -1, 0) + mut.list$Gain;

    col.data <- structure(
        c('firebrick1', 'white', 'dodgerblue'),
        names = c('Gain', 'Loss', 'SSM')
        );

    mutation.types <- get.mutation.type(colnames(mutation.data), clonal = FALSE);
    mut.labels <- get.mutation.labels(colnames(mutation.data));

    }

create.CNA.heatmap <- function(
    CNA.data,
    colours = c('dodgerblue', 'white', 'firebrick1'),
    same.as.matrix = TRUE,
    xaxis.tck = c(0, 0),
    yaxis.tck = c(.5, 0),
    ...
    ) {

    CNA.heatmap <- create.heatmap(
        x = CNA.data,
        clustering.method = 'none',
        print.colour.key = FALSE,
        colour.scheme = colours,
        at = c(-1, -0.5, 0, 1),
        total.colours = 4,
        axes.lwd = 1,
        same.as.matrix = same.as.matrix,
        xaxis.tck = xaxis.tck,
        yaxis.tck = yaxis.tck,
        ...
        );

    CNA.legend <- list(
        colours = colours[colours != 'white'],
        labels = c('Gain','Loss'),
        title = expression(bold(underline('Mutation Type')))
        );
    return(list(heatmap = CNA.heatmap, legend = CNA.legend));
    }