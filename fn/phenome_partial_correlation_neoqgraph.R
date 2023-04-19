require(qgraph)
require(igraph)
require(data.table)
require(cli)
require(gridExtra)
require(scales)

phenome_partial_correlation_neoqgraph <- function(
    x,
    from_var      = "phe1",
    to_var        = "phe2",
    cor_var       = "estimate",
    phe_group     = "neoplasms",
    savefile      = NULL,
    out_width     = 20,
    out_height    = 20,
    no_leafs      = FALSE,
    rollup_only   = FALSE,
    lauren_rollup = FALSE,
    qcut          = 0,
    corr_thresh   = NULL,
    plot_key      = FALSE
) {
    # prepare neoplasm phecodes
    pheinfo <- fread("https://gitlab.com/maxsal/public_data/-/raw/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
                    colClasses = "character",
                    showProgress = FALSE)
    if (no_leafs) { pheinfo <- pheinfo[leaf == 0, ] }
    if (rollup_only) { pheinfo <- pheinfo[rollup == 1, ] }
    if (lauren_rollup) {
        rollup    <- fread("https://gitlab.com/maxsal/public_data/-/raw/main/phewas/rollup.txt",
                           colClasses = "character",
                           showProgress = FALSE)[[1]]
        rollup_tf <- pheinfo[, phecode] %in% rollup
        pheinfo   <- pheinfo[rollup_tf, ]
    }
    if (!phe_group %in% pheinfo[, unique(group)]) { stop(paste0("'phe_group' (", phe_group, ") not in pheinfo groups: ",
                                                                paste0(pheinfo[, unique(group)], collapse = ", ")))}
    
    group_phe <- paste0("X", pheinfo[group == phe_group, phecode])
    if (phe_group == "neoplasms") {
        benign <- pheinfo[which(group == "neoplasms" &
            (grepl("[Bb]enign", description) |
             description %in% c(
                'Kaposi\'s sarcoma',
                'Lipoma',
                'Lipoma of skin and subcutaneous tissue',
                'Nevus, non-neoplastic',
                'Vascular hamartomas and non-neoplastic nevi',
                'Uterine leiomyoma',
                'Hemangioma and lymphangioma, any site',
                'Hemangioma of skin and subcutaneous tissue'))),
            phecode]
        treatment <- pheinfo[description %in% c('Bone marrow or stem cell transplant',
                                                'Radiotherapy', 'Chemotherapy',
                                                'Acquired absence of breast',
                                                'Screening for malignant neoplasms of the skin'),
                             phecode]
    }
    
    # prepare partial correlation data
    vars            <- c(from_var, to_var, cor_var)
    pcor_sub        <- as.data.table(x)[, ..vars]
    names(pcor_sub) <- c("from", "to", "Freq")
    if (!is.null(corr_thresh)) {
        pcor_sub <- pcor_sub[Freq >= corr_thresh, ]
    }
    pcor_sub        <- as.matrix(as_adj(graph.data.frame(pcor_sub[to %in% group_phe & from %in% group_phe, ]),
                                 type = "both", attr = "Freq"))
    colors <- rep('cadetblue3', ncol(pcor_sub))
    if (phe_group == "neoplasms") {
        colors[colnames(pcor_sub) %in% paste0("X", benign)]    <- "darkseagreen3"
        colors[colnames(pcor_sub) %in% paste0("X", treatment)] <- "gold2"
    }

    colnames(pcor_sub) <- gsub("X", "", colnames(pcor_sub))
    rownames(pcor_sub) <- gsub("X", "", rownames(pcor_sub))

    if (!is.null(savefile)) { pdf(file = paste0(gsub(".pdf", "", savefile), ".pdf"), width = out_width, height = out_height) }
    Q <- suppressWarnings(qgraph(as.data.frame(pcor_sub), layout = 'circle', shape = 'circle', legend = FALSE,
            label.cex = 1, node.width = 0.7,node.height = 0.7, 
            curveAll = TRUE, curveScale = FALSE, arrows= FALSE, 
            minimum = 0, maximum = 1, color = colors,
            negCol = 'red', posCol = 'gray42', cut = qcut))
    if  (!is.null(savefile)) {
        cli_alert_info("qgraph saved to {.path {savefile}}")
        dev.off()
    }
    
    # key
    if (phe_group == "neoplasms" & plot_key) {
        tt <- gridExtra::ttheme_default(core=list(bg_params = list(fill = 'gray95', col='black'),
                                               fg_params=list(fontface=3, fontsize=12)),
                                     colhead=list(fg_params=list(col="white", fontface=4L, fontsize=12), 
                                                  bg_params = list(fill = 'white', col='white')),
                                     rowhead=list(fg_params=list(col="white", fontface=4L, fontsize=12),
                                                  bg_params = list(fill = 'white', col='white')))

        find_cell <- function(table, row, col, name="core-fg"){
          l <- table$layout
          which(l$t==row & l$l==col & l$name==name)
        }
        DAT = data.frame(A = 'Malignant', B = 'Benign', C = 'Treatment-Related')
        SHORT = c('cadetblue3',  'darkseagreen3', 'gold2')
        ss <- tableGrob(DAT, rows = rep("", 1), theme = tt)
        for(i in 1:3){
          ind2 <- find_cell(ss, 2,i+1, "core-bg")
          ss$grobs[ind2][[1]][["gp"]] <- grid::gpar(fill=SHORT[i], col= 'black')
        }

        if (!is.null(savefile)) { pdf(file = paste0(gsub(".pdf", "", savefile), "_key.pdf"), width = 6, height = 3) }
        grid.arrange(ss, ncol = 1)
        if (!is.null(savefile)) { dev.off() }
    }
    
}
