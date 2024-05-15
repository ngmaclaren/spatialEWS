## Directed networks: yeast, air_carrier (maybe), jung
## Weighted: product export,
## Main euroroad
## -------------
## Konect: iceland, Wikipedia talk (ht), hamsterster (I think the lcc is smaller), contiguous USA, gene fusion

library(igraph)
set.seed(123)

                                        # Only needed once...
## usair <- read.csv("./us_air_traffic/edges.csv")
## usair <- usair[usair$year == 2020, ]
## write.csv(usair, file = "./us_air_traffic/edges.csv", row.names = FALSE)

dirs <- c(
    "adjnoun",
    "canton",
    "catlins",
    "chesapeake",
    "contigusa",
    "dolphin",
    "drug",
    "ecoli",
    "email_company",
    "euroroad",
    "flamingo",
    "football",
    "gap_junction_herm",
    "genefusion",
    "hk100",
    "iceland",
    "jazz",
    "jung-c",
    "metabolic",
    "montreal",
    "netsci",
    "physician_trust",
    "protein",
    "proximity",
    "SITC",
    "students",
    "train_terrorists",
    "us_air_traffic",
    "wiki-ht",
    "windsurfers",
    "yeast"
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
## networks$tree <- make_tree(N, mode = "undirected")
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
