library(BoutrosLab.plotting.general);
#x = data.frame with e.g. genes in rows, and feature of interest in columns
#should be a data.frame of logicals, where TRUE = some colour and FALSE = white, NA = grey
#colour.scheme = names list of colours, where names match columns in data.frame, if not provided, black will be used for all columns
#cluster.rows
#FALSE.colour = 'white'
create.discrete.colour.heatmap <- function(
    x,
    colour.scheme,
    cluster.rows = TRUE,
    FALSE.colour = 'white',
    NA.colour = 'darkgray',
    key.title = NULL,
    ...
    ) {

    # If more than one term, perform clustering
    if (ncol(x) > 1 & cluster.rows) {
        x.order <- hclust(dist(t(x)))$order;
        x <- x[, x.order];
        colour.scheme <- colour.scheme[x.order];
        #convert to colour
        x.colours <- t(data.frame(sapply(names(x), function(column)
            ifelse(x[, column], colour.scheme[[column]], FALSE.colour)
            )));
        } else {
            x.colours <- data.frame(sapply(names(x), function(column)
                ifelse(x[, column], colour.scheme[[column]], FALSE.colour)
                ));
        }
    if (any(is.na(x.colours))) {
        x.colours[is.na(x.colours)] <- NA.colour;
        }
    #create legend
    if (!is.expression(key.title)) {
        key.title <- bquote(bold(underline(.(key.title))));
        }
    x.legend <- list(
        colours = as.character(colour.scheme),
        labels = names(colour.scheme),
        title = key.title,
        cex = .6
        );

    x.heatmap <- create.heatmap(
        x = x.colours,
        clustering.method = 'none',
        cluster.dimensions = 'none',
        print.colour.key = FALSE,
        input.colours = TRUE,
        yaxis.tck = 0,
        xaxis.tck = 0,
        axes.lwd = 1,
        ...
        );
    #return heatmap and legend
    return(list(heatmap = x.heatmap, legend = x.legend));
    }