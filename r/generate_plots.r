# libraries and such ------------------------------------------------------
library(ms)
library(qs)
library(data.table)
library(glue)
library(patchwork)

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

prevs <- list() 
for (i in seq_along(bbs)) {
    prevs[[i]] <- qread(glue("{data_path}{bbs[i]}_prevs.qs"))
    names(prevs)[i] <- bbs[i]
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
