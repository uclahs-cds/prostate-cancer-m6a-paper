#diverging colour scheme for data centered around 0
create.diverging.heatmap <- function(
    data,
    add.cell.text = FALSE,
    key.min = -3,
    key.max = 3,
    mid.point = 0,
    key.colour.interval.num = 50,
    colours = c('darkorchid4', 'white', 'darkgreen'),
    filename = NULL,
    same.as.matrix = FALSE,
    scale.data = 'none',
    key.title = '',
    xaxis.tck = c(0, 0),
    yaxis.tck = c(.3, 0),
    return.plot = TRUE,
    ...
    ) {

    if ('row' == scale.data) {
        data <- t(scale(t(data)));
        } else if ('col' == scale.data) {
        data <- scale(data);
        }

    key.scale <- c(
        seq(key.min, mid.point, -key.min / key.colour.interval.num),
        seq(mid.point, key.max, key.max / key.colour.interval.num)
        );
    key.scale <- unique(key.scale);
    if (add.cell.text & same.as.matrix) {
        #stop('add.cell.text and same.as.matrix not currently implemented together');
        cell.text = round(data, digits = 2);
        row.pos = rep(rev(1:nrow(data)), ncol(data));
        col.pos = as.vector(sapply(1:ncol(data), function(row) rep(row, nrow(data))));
        } else if (add.cell.text) {
        cell.text = round(data, digits = 2);
        col.pos = rep(1:nrow(data), ncol(data));
        row.pos = as.vector(sapply(1:ncol(data), function(row) rep(row, nrow(data))));
        } else {
        cell.text = '';
        row.pos = NULL;
        col.pos = NULL;
        }
    heatmap <- BoutrosLab.plotting.general::create.heatmap(
        x = data,
        filename = filename,
        at = key.scale,
        colour.scheme = colours,
        print.colour.key = FALSE,
        scale.data = FALSE,
        cell.text = cell.text,
        row.pos = row.pos,
        col.pos = col.pos,
        same.as.matrix = same.as.matrix,
        xaxis.tck = xaxis.tck,
        yaxis.tck = yaxis.tck,
        axes.lwd = 1,
        ...
        );
    if (!is.expression(key.title)) {
        key.title <- expression(bold(underline('Z'[key.title])));
        }
    if (return.plot) {
        heatmap.legend <- list(
            colours = colours,
            labels = c(
                as.expression(bquote(''<=.(round(key.scale[1], digits = 2)))),
                as.character(median(key.scale)),
                as.expression(bquote(''>=.(round(key.scale[length(key.scale)], digits = 2))))
                ),
            title = key.title,
            continuous = TRUE,
            tck = 0,
            at = c(10, 50, 90),
            cex = .6
            );
        return(list(heatmap = heatmap, legend = heatmap.legend));
        }
    }
