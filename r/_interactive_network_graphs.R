library(ms)
library(qs)
library(data.table)
library(tidyverse)
library(igraph)
library(networkD3)
library(geomnet)
library(visNetwork)

# Load data
aou <- qread("data/aou_phenome_partial_correlations_20230309.qs")
aou <- aou[abs(estimate) > 0.3, ]

# networkD3 example
# simpleNetwork(
#     aou[, .(from = phe1, to = phe2)],
#     height = "500px", width = "800px",
#     Source = 1,
#     Target = 2,
#     linkDistance = 10,
#     charge = -900,
#     fontSize = 14,
#     fontFamily = "sans-serif",
#     linkColour = "#ccc",
#     nodeColour = "#69b3a2",
#     opacity = 0.9,
#     zoom = TRUE
# )

nodes <- data.table(
    label = unique(c(aou[, phe1], aou[, phe2]))
)[order(as.numeric(gsub("X", "", label))), id := 0:(.N - 1)][]
edges <- aou[, .(from_label = phe1, to_label = phe2, width = estimate)]
edges <- merge.data.table(
    edges,
    nodes[, .(from_id = id, label)],
    by.x = "from_label",
    by.y = "label"
)
edges <- merge.data.table(
    edges,
    nodes[, .(to_id = id, label)],
    by.x = "to_label",
    by.y = "label"
)

cols <- merge.data.table(
    nodes[, .(label = gsub("X", "", label), id)],
    ms::pheinfo[, .(phecode, color)],
    by.x = "label",
    by.y = "phecode"
    )


(net <- forceNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "from_id",
    Target = "to_id",
    Value = "width",
    NodeID = "label",
    Group = "id",
    opacity = 0.9
    )
)

tmp <- data.table(
    from = edges[, from_label],
    to = edges[, to_label],
    Freq = 1,
    to_description = merge.data.table(
        edges[, .(to_label)],
        ms::pheinfo[, .(phecode = paste0("X", phecode), description)],
        by.x = "to_label", by.y = "phecode")[, description],
    from_description = merge.data.table(
        edges[, .(from_label)],
        ms::pheinfo[, .(phecode = paste0("X", phecode), description)],
        by.x = "from_label", by.y = "phecode"
    )[, description],
    to_category = merge.data.table(
        edges[, .(to_label)],
        ms::pheinfo[, .(phecode = paste0("X", phecode), group)],
        by.x = "to_label", by.y = "phecode"
    )[, group]
)

merge.data.table(
    data.table(phecode = gsub("X", "", ORDERING)),
    ms::pheinfo[, .(phecode, group)]
)[, group]

net = network(tmp, directect = FALSE, matrix.type = "edgelist")
ORDERING = net %v% "vertex.names"
net %v% "Category" = merge.data.table(
    data.table(phecode = gsub("X", "", ORDERING)),
    ms::pheinfo[, .(phecode, group)]
)[, group]
net %v% "Description" = as.character(c(unique(tmp[, to_description]), unique(tmp[, from_description])))
net %v% "Size" = 1
net %v% "CodeName" = c(rep("", length.out = length(unique(tmp[, to])), substring(as.character(unique(tmp[, from])), 2)))

cols <- ms::pheinfo[, .(group, color)]
cols_pal <- cols[, color]
names(cols_pal) <- cols[, group]

(net_plot <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend, label = Description)) +
    geom_edges(color = "grey50") +
    geom_nodes(mapping = aes(color = Category, size = Size)) +
    theme_blank() +
    guides(size = FALSE) +
    geom_nodetext(aes(label = CodeName, fontface = "bold", size = 1)) +
    scale_color_manual(values = cols_pal) +
    # scale_color_manual(values = c("gold2", "purple", "cadetblue3")) +
    scale_size(range = c(1, 5)))

(net_plotly <- ggplotly(net_plot, tooltip = c("vertex.names", "Description", "Category"),
dynamicTicks = FALSE, layerDate = 1, originalData = TRUE, height = 800, width = 1200))

htmlwidgets::saveWidget(net_plotly, "~/Downloads/test_net.html")

# visNetwork example
nodes <- as.data.frame(data.table(
    label = unique(c(aou[, phe1], aou[, phe2]))
)[order(as.numeric(gsub("X", "", label))), id := 0:(.N - 1)][])

edges <- as.data.frame(aou[, .(from = phe1, to = phe2, width = estimate)])

graph <- graph_from_data_frame(edges, directed = FALSE)
cluster <- cluster_louvain(graph)
cluster_df <- data.frame(as.list(membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- rownames(cluster_df)

nodes <- left_join(nodes, cluster_df, by = "label")
colnames(nodes)[3] <- "group"

visNetwork(nodes, edges)
