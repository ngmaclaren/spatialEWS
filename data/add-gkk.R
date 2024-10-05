library(igraph)
set.seed(123)

N <- 100
M <- 300

g <- largest_component(sample_fitness_pl(N, M, 2))

vcount(g)
ecount(g)
plot(g, vertex.size = 5, vertex.label = "")

networks <- readRDS("networks-pregkk.rds")
networks <- c(networks, list(gkk = g))
saveRDS(networks, "networks.rds")
