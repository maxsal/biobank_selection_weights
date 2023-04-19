require(PheWAS)
require(ComplexHeatmap)
require(igraph)
require(cli)

phenome_partial_correlation_heatmap <- function(
    x,
    savefile   = NULL,
    out_width  = 20,
    out_height = 20,
    plot_title = NULL,
    attr_col   = "estimate"
) {
    if (!attr_col %in% names(x)) {stop(paste0("'attr_col' (", attr_col, ") not in names of 'x'. change 'attr_col'"))}
    mat <- as.matrix(get.adjacency(graph.data.frame(as.data.frame(x), directed = FALSE), type = "both", attr = attr_col))
    if (is.null(plot_title)) { plot_title <- "Correlations" }

    color   <- rep(NA, length(colnames(mat)))
    groups  <- rep(NA, length(colnames(mat)))
    phecode <- rep(NA, length(colnames(mat)))
    gl      <- list()
    pheinfo <- PheWAS::pheinfo
    for (i in 1:length(colnames(mat))) {
        color[i]   <- pheinfo$color[paste0("X", pheinfo$phecode) == colnames(mat)[i]]
        groups[i]  <- pheinfo$group[paste0("X", pheinfo$phecode) == colnames(mat)[i]]
        phecode[i] <- pheinfo$phecode[paste0("X", pheinfo$phecode) == colnames(mat)[i]] 
        gl[[i]]    <- paste0("X", phecode[i])
    }
    names(gl)          <- groups
    phe_col            <- unique(as.data.table(pheinfo)[, .(color, group)])
    color_named        <- phe_col[, color]
    names(color_named) <- phe_col[, group] 

    gd <- structure(rep(names(gl), times = sapply(gl, length)), names = unlist(gl))
    gd <- gd[colnames(mat)]
    
    if (!is.null(savefile)) { pdf(file = savefile, width = out_width, height = out_height) }
    col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
    ht <- Heatmap(mat, name = plot_title,  #top_annotation = ha,
            col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, 
            show_row_names = FALSE, show_column_names = FALSE, 
            top_annotation = HeatmapAnnotation(group = gd, col = list(group = color_named), show_legend = FALSE),
            right_annotation = HeatmapAnnotation(group = gd, col = list(group = color_named), which = 'row'))
    draw(ht)
    if (!is.null(savefile)) {
        cli_alert_info("heatmap saved to {.path {savefile}}")
        dev.off()
    }
}
