library(igraph)

dirs <- c(
    "adjnoun",
    "canton",
    "catlins",
    "chesapeake",
    "dolphin",
    "drug",
    "ecoli",
    "email_company",
    "flamingo",
    "gap_junction_herm",
    "jazz",
    "metabolic",
    "montreal",
    "netsci",
    ## "pdzbase",
    "physician_trust",
    "protein",
    "proximity",
    "students",
    "train_terrorists",
    "windsurfers"
)

gather_network <- function(dir) {
    el <- read.csv(paste0("./", dir, "/edges.csv"))
    g <- simplify(largest_component(graph_from_data_frame(el, directed = FALSE)))
    if("weight" %in% edge_attr_names(g)) g <- delete_edge_attr(g, "weight")
    g$name <- dir
    return(g)
}

networks <- lapply(dirs, gather_network)
names(networks) <- dirs

N <- 100
networks$smallworld <- sample_smallworld(1, N, 4, 0.02)
networks$barabasialbert <- sample_pa(N, m = 2, directed = FALSE, start.graph = make_full_graph(3))
networks$erdosrenyi <- largest_component(sample_gnp(N, 0.05))
networks$tree <- make_tree(N, mode = "undirected")
networks$lattice <- make_lattice(length = 10, dim = 2, circular = TRUE)

networks <- networks[order(sapply(networks, vcount))]

sapply(networks, vcount)
sapply(networks, ecount)
sapply(networks, is_simple) # students network is not simple. 
sapply(networks, is_connected)
sapply(networks, is_directed)
sapply(networks, is_weighted)

## g <- networks$students
## plot(g, vertex.size = 5, vertex.label = "")
.plot <- function(g) plot(g, vertex.size = 5, vertex.label = "", main = g$name)
for(g in networks) {
    dev.new()
    .plot(g)
}

saveRDS(networks, "./networks.rds")




## Code to import catlins
## catlins <- largest_component(
##     ## graph_from_adjacency_matrix(
##     graph_from_biadjacency_matrix(
##         as.matrix(read.table("./catlins/catlins-matrix.txt", sep = "\t")), directed = FALSE
##     )
## )
## df <- as_data_frame(catlins, "edges")
## colnames(df) <- c("source", "target")
## write.csv(df, file = "./catlins/edges.csv", row.names = FALSE)
