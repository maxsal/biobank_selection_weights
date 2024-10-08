suppressPackageStartupMessages({
  library(ggnetwork)
  library(network)
  library(data.table)
  library(cli)
})

phenome_partial_correlation_network <- function(
    x,
    from_var   = "var1",
    to_var     = "var2",
    cor_var    = "cor",
    prev_var   = "prev_unweighted",
    thresh     = 0.3,
    savefile   = NULL,
    out_width  = 12,
    out_height = 8,
    plot_title = NULL,
    prevs,
    show_color_legend = NULL,
    show_size_legend  = NULL
) {
    # initialize
    vars         <- c(from_var, to_var, cor_var)
    edges        <- as.data.table(x)[, ..vars]
    names(edges) <- c("from", "to", "Freq")
    edges        <- as.data.frame(edges[!is.na(Freq), ][abs(Freq) >= thresh, ])
    cli_alert_info(paste0("plotting correlations with absolute value >= ", thresh, " (n = ", nrow(edges), ")"))
    edges[, 1] <- as.character(edges[, 1])
    edges[, 2] <- as.character(edges[, 2])
    if (is.null(plot_title)) plot_title <- "Correlations"

    net <- network(edges, directed = FALSE, matrix.type = 'edgelist')
    network::set.edge.attribute(net, "Correlation", edges$Freq)
    nodecol = groupcol = prevcol = cormag = vertexnames = codenames = c()
    phecode_info <- as.data.table(ms::pheinfox)
    cli_progress_bar(name = "preparing data for network plot", total = length(net$val))
    for(i in seq_along(net$val)){
        code        <- net$val[[i]]$vertex.names
        nodecol     <- c(nodecol, phecode_info[phecode_info$phecode == code, ]$color)
        groupcol    <- c(groupcol, phecode_info[phecode_info$phecode == code, ]$group)
        prevcol     <- c(prevcol, prevs[phecode == code, get(prev_var)])
        vertexnames <- c(vertexnames, phecode_info$description[phecode_info$phecode == code])
        codenames   <- c(codenames, code)
        cli_progress_update()
    }
    cli_progress_done()

    net %v% "Category"    <- groupcol
    net %v% "colors"      <- nodecol
    net %v% "Prevalence"  <- round(prevcol,4)
    net %v% "Description" <- paste0(vertexnames, ' (', codenames, ')')
    net %v% "Phecode"     <- codenames

    p <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend, label = Description)) +
        geom_edges(color = "grey50" ) +
        #geom_edges(aes(size =edgeval/100), show.legend= FALSE) +
        geom_nodes(mapping = aes(color = Category, size =  Prevalence)) +
        scale_color_manual(values = unique(nodecol)) +
        guides(
            col  = show_color_legend,
            size = show_size_legend
        ) +
        labs(
            color   = "Disease Category",
            title   = plot_title,
            caption = paste0(format(nrow(edges), big.mark = ","), " across ", format(length(unique(c(edges$to, edges$from))), big.mark = ","),
                             " unique phecodes with correlations greater than ", thresh)
        ) +
        theme_blank() +
        theme(
                legend.text = element_text(size = 8, colour = "black"),
                legend.key.size = unit(0.001, 'cm'),
                plot.caption = element_text(hjust = 0)
            )
    
    if (!is.null(savefile)) {
        ggsave(
            plot     = p,
            filename = savefile,
            width    = out_width,
            height   = out_height
        )
        cli_alert_success("plot saved to {.path {savefile}}")
    }

    return(p)

}
