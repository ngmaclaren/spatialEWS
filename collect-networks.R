library(igraph)
                                        # for interactive use
.plot <- function(g, ...) plot(g, vertex.size = 4, vertex.label = "", ...)

## Naming convention: initials/abbreviations for model networks, written-out names for empirical networks
## Currently have 42. Would like a rounder number, like 40 or 45. 

#### MODEL NETWORKS ####
N <- 100

## deterministic
pl <- make_lattice(length = floor(sqrt(N)), dim = 2, circular = TRUE)
sl <- make_lattice(length = floor(sqrt(N)), dim = 2, circular = FALSE)

t2 <- make_tree(N, 2, "undirected")
t3 <- make_tree(N, 3, "undirected")
t4 <- make_tree(N, 4, "undirected")

## stochastic
er <- largest_component(sample_gnm(N, 2*N))

p <- 0.05
sw2 <- sample_smallworld(dim = 1, size = N, nei = 2, p = p)
sw3 <- sample_smallworld(dim = 1, size = N, nei = 3, p = p)
sw4 <- sample_smallworld(dim = 1, size = N, nei = 4, p = p)

ba2 <- with(list(m = 2), sample_pa(N, m = m, directed = FALSE, start.graph = make_full_graph(m + 1)))
ba3 <- with(list(m = 3), sample_pa(N, m = m, directed = FALSE, start.graph = make_full_graph(m + 1)))
ba4 <- with(list(m = 4), sample_pa(N, m = m, directed = FALSE, start.graph = make_full_graph(m + 1)))

#### EMPIRICAL NETWORKS ####
            
### empirical
dolphin <- graph_from_edgelist(as.matrix(read.csv("../networks/dolphin.csv")), FALSE)
proximity <- graph_from_edgelist(as.matrix(read.csv("../networks/proximity.csv")), FALSE)
genefusion <- graph_from_edgelist(as.matrix(read.csv("../networks/gene-fusion.csv")), FALSE)
jazz <- graph_from_edgelist(as.matrix(read.csv("../networks/jazz.csv")), FALSE)
netsci <- graph_from_edgelist(as.matrix(read.csv("../networks/netsci.csv")), FALSE)
pdzbase <- graph_from_edgelist(as.matrix(read.csv("../networks/pdzbase.csv")), FALSE)
## ecological networks
chesapeake <- graph_from_edgelist(as.matrix(read.csv("../networks/chesapeake.csv")), FALSE)
joern <- graph_from_biadjacency_matrix(read.table("../networks/originals/joern_1979_altuda.txt", sep = "\t"))
parsnip <- graph_from_biadjacency_matrix(read.table("../networks/originals/parsnip_p.txt", sep = "\t"))
huron <- graph_from_biadjacency_matrix(read.table("../networks/originals/sbay_huron_p.txt", sep = "\t"))
akatoreA <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/AkatoreAtxt.txt", sep = "\t"))
akatoreB <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/AkatoreBtxt.txt", sep = "\t"))
berwick  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Berwicktxt.txt", sep = "\t"))
blackrock  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Blackrocktxt.txt", sep = "\t"))
broad  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Broadtxt.txt", sep = "\t"))
catlins  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Catlinstxt.txt", sep = "\t"))
coweeta17  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Coweeta17txt.txt", sep = "\t"))
coweeta1  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Coweeta1txt.txt", sep = "\t"))
dempstersA  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/DempstersAutxt.txt", sep = "\t"))
dempstersS  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/DempstersSptxt.txt", sep = "\t"))
dempstersSu  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/DempstersSutxt.txt", sep = "\t"))
german  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Germantxt.txt", sep = "\t"))
healy  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Healytxt.txt", sep = "\t"))
kyeburn  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Kyeburntxt.txt", sep = "\t"))
lilkyeburn  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/LilKyeburntxt.txt", sep = "\t"))
martins  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Martinstxt.txt", sep = "\t"))
narrowdale  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Narrowdaletxt.txt", sep = "\t"))
northcol  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/NorthColtxt.txt", sep = "\t"))
stony  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Stonytxt.txt", sep = "\t"))
suttonAu  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/SuttonAutxt.txt", sep = "\t"))
suttonSp  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/SuttonSptxt.txt", sep = "\t"))
suttonSu  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/SuttonSutxt.txt", sep = "\t"))
troy  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Troytxt.txt", sep = "\t"))
venlaw  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Venlawtxt.txt", sep = "\t"))
## canton  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Cantontxt.txt", sep = "\t"))
## powder  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Powdertxt.txt", sep = "\t"))

                                        # Fragile, but works if session is empty/clear
networknames <- ls()[-which(ls() %in% c("N", "p"))]
gl <- lapply(networknames, get) 
names(gl) <- networknames
                                        # then loop through, replacing each with simplify(largest_component(g))
for(i in seq_along(gl)) {
    gl[[i]] <- simplify(largest_component(gl[[i]]))
}

## for checking
## for(i in seq_along(gl)) {dev.new(); .plot(gl[[i]], main = names(gl[i]))}

saveRDS(gl, file = "./data/networks.rds")
