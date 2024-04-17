library(igraph)

dirs <- c("adjnoun", "dolphin", "montreal", "students")

gather_network <- function(dir) {
    el <- read.csv(paste0("./", dir, "/edges.csv"))
    g <- simplify(largest_component(graph_from_data_frame(el, directed = FALSE)))
}

networks <- lapply(dirs, gather_network)
names(networks) <- dirs

N <- 100
networks$smallworld <- sample_smallworld(1, N, 4, 0.01)
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
saveRDS(networks, "./networks.rds")
