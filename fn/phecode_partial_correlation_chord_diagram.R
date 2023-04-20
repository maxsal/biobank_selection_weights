suppressPackageStartupMessages({
  library(GGally)
  library(network)
  library(ComplexHeatmap)
  library(circlize)
})

phenome_partial_correlation_chord_diagram <- function(
    x,
    from_var   = "phe1",
    to_var     = "phe2",
    cor_var    = "estimate",
    thresh     = 0.2,
    savefile   = NULL,
    out_width  = 10,
    out_height = 10,
    plot_title = NULL
    ) {
    # initialize
    vars <- c(from_var, to_var, cor_var)
    edges <- as.data.table(x)[, ..vars]
    names(edges) <- c("from", "to", "Freq")
    edges <- as.data.frame(edges[!is.na(Freq), ][abs(Freq) >= thresh, ])
    message(paste0("plotting correlations with absolute value >= ", thresh, " (n = ", nrow(edges), ")"))
    edges[, 1] = as.character(edges[, 1])
    edges[, 2] = as.character(edges[, 2])

    if (is.null(plot_title)) { plot_title <- "Correlations" }
    
    # prep phecode info
    phecode_info  <- as.data.table(PheWAS::pheinfo)
    GROUPS = unique(phecode_info$group)
    COLOR = rep(NA,length(GROUPS))
    for(i in 1:length(GROUPS)){
    COLOR[i] = unique(phecode_info$color[phecode_info$group == GROUPS[i]])
    }
    COLOR[COLOR == 'black'] = 'purple'
    COLOR[COLOR == 'darkblue'] = 'cornflowerblue'
    COLOR[COLOR == 'gray44'] = 'azure2'
    
    # rename phecodes for plotting
    for(i in 1:length(GROUPS)){
        code <- paste0('X', phecode_info[group == GROUPS[i], phecode])
        edges$from[edges$from %in% code] = GROUPS[i]
        edges$to[edges$to %in% code] = GROUPS[i]
    }
    edges = edges[edges$from != edges$to,]

    col_fun = colorRamp2(c(-1,-0.5, 0, 0.5,1), c('red',"red", "white", "black", 'black'), transparency = c(0,0,1,0,0))
    weights = rep(NA,length(edges$Freq))
    weights[edges$Freq < 0.4] = 1
    weights[edges$Freq >= 0.4 & edges$Freq < 0.5] = 5
    weights[edges$Freq >= 0.5 & edges$Freq < 0.6] = 10
    weights[edges$Freq >= 0.6 & edges$Freq < 0.7] = 50
    weights[edges$Freq >= 0.7 & edges$Freq < 0.8] = 100
    weights[edges$Freq >= 0.8 & edges$Freq < 0.9] = 500
    weights[edges$Freq >= 0.9 & edges$Freq <= 1] = 1000

    # plot
    if (!is.null(savefile)) { pdf(file = savefile, width = out_width, height = out_height) }
    circos.par(gap.degree = as.numeric(1/length(edges$Freq)))
    chordDiagram(edges, grid.col = COLOR[which(GROUPS %in% c(edges$to, edges$from))], col = col_fun(edges$Freq),
                annotationTrack = c('grid'), 
                annotationTrackHeight = convert_height(c(30), "mm"), scale = FALSE,
                link.lwd = weights, link.largest.ontop = TRUE
    )
    lgd_links = ComplexHeatmap::Legend(at = c(-1, 0, 1),col_fun = col_fun, title_position = "topleft",
    title = plot_title)
    ComplexHeatmap::draw(lgd_links, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
    for(si in get.all.sector.index()) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
                facing = "clockwise", niceFacing = TRUE, col = 'black', cex = 1,adj = 0.5)
    }
    circos.clear()
    if (!is.null(savefile)) {
        dev.off()
        message(paste0("chord diagram saved to ", savefile))
        }
}