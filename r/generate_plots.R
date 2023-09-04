# libraries and such ------------------------------------------------------
ms::libri(
    ms, qs, data.table, glue, patchwork, ggplot2
)

source("~/Dropbox (University of Michigan)/projects/dissertation/aim_one/fn/phenome_partial_correlation_network.R")
data_path <- "~/Dropbox (University of Michigan)/projects/dissertation/all_of_us/data/"

# load data ---------------------------------------------------------------
bbs      <- c("mgi", "ukb", "aou")
bb_dates <- c("20220822", "20221117", "20230309")

parcors <- list()
for (i in seq_along(bbs)) {
    parcors[[i]] <- qread(glue("{data_path}{bbs[i]}_phenome_partial_correlations_{bb_dates[i]}.qs"))
    names(parcors)[i] <- bbs[i]
}

#setnames(parcors[["mgi"]], c("var1", "var2", "cor"), c("phe1", "phe2", "estimate"))
#setnames(parcors[["ukb"]], c("var1", "var2", "cor"), c("phe1", "phe2", "estimate"))
#setnames(parcors[["aou"]], c("var1", "var2", "cor"), c("phe1", "phe2", "estimate"))

wparcors <- list()
for (i in seq_along(bbs[bbs])) {
    wparcors[[i]] <- qread(glue("{data_path}{bbs[i]}_phenome_weighted_partial_correlations_{bb_dates[i]}.qs"))
    names(wparcors)[i] <- bbs[i]
}
#setnames(wparcors[["mgi"]], c("var1", "var2", "cor"), c("phe1", "phe2", "estimate"))
#setnames(wparcors[["ukb"]], c("var1", "var2", "cor"), c("phe1", "phe2", "estimate"))
#setnames(wparcors[["aou"]], c("var1", "var2", "cor"), c("phe1", "phe2", "estimate"))

prevs <- list() 
for (i in seq_along(bbs)) {
    prevs[[i]] <- qread(glue("{data_path}{bbs[i]}_prevs.qs"))
    names(prevs)[i] <- bbs[i]
}

common_phecodes <- lapply(prevs, \(x) x[n >= 20, phecode]) |>
    Reduce(intersect, x = _)
    
for (i in seq_along(bbs)) {
    parcors[[i]] <- parcors[[i]][var1 %in% common_phecodes & var2 %in% common_phecodes]
    wparcors[[i]] <- wparcors[[i]][var1 %in% common_phecodes & var2 %in% common_phecodes]
}

# partial correlation graphs ----------------------------------------------
plots <- list()
for (i in seq_along(bbs)) {
    plots[[i]] <- phenome_partial_correlation_network(
        x = parcors[[i]], prevs = prevs[[i]], savefile = glue("results/{bbs[i]}_network.pdf"),
        plot_title = glue("Correlation network in {toupper(bbs[i])}"),
        show_color_legend = "none"
    )
    names(plots)[i] <- bbs[i]
}

col_legend <- cowplot::get_legend(
        phenome_partial_correlation_network(
        x = parcors[[1]], prevs = prevs[[1]],
        plot_title = glue("Correlation network in {toupper(bbs[i])}"),
        show_color_legend = NULL,
        show_size_legend = "none"
    ) +
    theme(legend.position = "bottom")
)

wplots <- list()
for (i in seq_along(bbs)) {
    wplots[[i]] <- phenome_partial_correlation_network(
        x = wparcors[[i]], prevs = prevs[[i]], savefile = glue("results/{bbs[i]}_wnetwork.pdf"),
        prev_var = "prev_weighted",
        plot_title = glue("Correlation network in {toupper(bbs[i])}"),
        show_color_legend = "none"
    )
    names(wplots)[i] <- bbs[i]
}

# stack plots -------------------------------------------------------------
stacked <-  ((((plots[["aou"]] +
    labs(
        title = "A. All of Us"
    ) + theme(
        plot.title = element_text(face = "bold")
    )) /
    (plots[["mgi"]] +
        labs(
            title = "B. Michigan Genomics Initiative"
        ) + theme(
            plot.title = element_text(face = "bold")
        )) /
    (plots[["ukb"]] +
        labs(
            title = "C. UK Biobank"
        ) + theme(
            plot.title = element_text(face = "bold")
        )))) /
        (col_legend)) +
        plot_layout(heights = c(3, 3, 3, 1))

ggsave(
    stacked,
    filename = "results/stacked_networks.pdf",
    width = 11.75, height = 15, device = cairo_pdf
)

wstacked <- ((((wplots[["aou"]] +
        labs(
            title = "A. All of Us"
        ) + theme(
            plot.title = element_text(face = "bold")
        )) /
        (wplots[["mgi"]] +
        labs(
            title = "B. Michigan Genomics Initiative"
        ) + theme(
            plot.title = element_text(face = "bold")
        )) /
    (wplots[["ukb"]] +
        labs(
            title = "C. UK Biobank"
        ) + theme(
            plot.title = element_text(face = "bold")
        )))) /
    (col_legend)) +
    plot_layout(heights = c(3, 3, 3, 1))

ggsave(
    wstacked,
    filename = "results/stacked_wnetworks.pdf",
    width = 11.75, height = 10, device = cairo_pdf
)

aou_networks <- ((plots[["aou"]] +
    labs(
        title = "A. Unweighted"
    ) +
    theme(plot.caption = element_blank())) / (wplots[["aou"]] +
    labs(
        title = "B. Weighted"
    ) +
    theme(plot.caption = element_blank())) /
    col_legend) +
    plot_annotation(
        title = "Partial correlation networks in All of Us"
    ) +
        plot_layout(heights = c(3, 3, 1)) &
    theme(plot.title = element_text(face = "bold"))

mgi_networks <- ((plots[["mgi"]] +
    labs(
        title = "A. Unweighted"
    ) +
    theme(plot.caption = element_blank())) / (wplots[["mgi"]] +
    labs(
        title = "B. Weighted"
    ) +
    theme(plot.caption = element_blank())) /
    col_legend) +
    plot_annotation(
        title = "Partial correlation networks in the Michigan Genomics Initiative"
    ) +
    plot_layout(heights = c(3, 3, 1)) &
    theme(plot.title = element_text(face = "bold"))

ukb_networks <- ((plots[["ukb"]] +
    labs(
        title = "A. Unweighted"
    ) +
    theme(plot.caption = element_blank())) / (wplots[["ukb"]] +
    labs(
        title = "B. Weighted"
    ) +
    theme(plot.caption = element_blank())) /
    col_legend) +
    plot_annotation(
        title = "Partial correlation networks in the UK Biobank"
    ) +
    plot_layout(heights = c(3, 3, 1)) &
    theme(plot.title = element_text(face = "bold"))

ggsave(
    aou_networks,
    filename = "results/aou_networks.pdf",
    width = 11.75, height = 10, device = cairo_pdf
)

ggsave(
    mgi_networks,
    filename = "results/mgi_networks.pdf",
    width = 11.75, height = 10, device = cairo_pdf
)

ggsave(
    ukb_networks,
    filename = "results/ukb_networks.pdf",
    width = 11.75, height = 10, device = cairo_pdf
)
