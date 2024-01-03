ms::libri(
    ms, data.table, ggplot2
)

d <- fread("results/hypertension_female_mgi20220822_ukb20221117.csv")
d[, `:=` (
    or = exp(beta),
    weight = fcase(
        is.na(weight) | weight == "", "Unweighted",
        grepl("ps", weight), "PostStrat",
        grepl("ip", weight) | weight == "weight", "IPW",
        default = NA
    )
    )]

d |>
    ggplot(aes(x = dataset, y = or)) +
    geom_hline(yintercept = 1, color = "black", linewidth = 1) +
    geom_point(aes(color = weight, shape = weight), size = 4) +
    ylim(0, NA) +
    labs(
        title = "Female odds ratio for hypertension",
        x = "Dataset",
        y = "Odds ratio"
    ) +
    scale_color_ms() +
    theme_ms()
