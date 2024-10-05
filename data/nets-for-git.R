library(igraph)

networks <- readRDS("networks.rds")

models <- c("lattice", "erdosrenyi", "smallworld", "barabasialbert", "hk100", "gkk")

modelnetworks <- networks[models]

## lapply(modelnetworks, function(g) {dev.new(); plot(g, vertex.size = 5, vertex.label = "")})

saveRDS(modelnetworks, "modelnetworks.rds")
