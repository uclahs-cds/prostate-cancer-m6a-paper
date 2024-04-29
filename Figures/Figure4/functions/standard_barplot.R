
create.standard.barplot <- function(
    data,
    x,
    xaxis.lab = rep('', nrow(buffa)),
    yaxis.tck = c(0.7, 0),
    xaxis.tck = c(0, 0),
    yaxis.cex = 0.8,
    xlab.label = '',
    ylab.cex = .7,
    col = 'black',
    box.ratio = .35,
    ylab.label = '',
    ...
    ){
    if (!is.expression(ylab.label)) {
        ylab.label <- bquote(bold(.(ylab.label)));
        }
    data$x <- data[, x];
    data$y <- 1:nrow(data);
    barplot <- create.barplot(
        data = data,
        formula = x ~ y,
        style = 'Nature',
        ylab.label = ylab.label,
        ylab.cex = ylab.cex,
        yaxis.cex = yaxis.cex,
        xlab.label = xlab.label,
        yaxis.tck = yaxis.tck,
        xaxis.tck = xaxis.tck,
        xaxis.lab = xaxis.lab,
        col = col,
        box.ratio = box.ratio,
        ...
        );

    return(barplot)
    }