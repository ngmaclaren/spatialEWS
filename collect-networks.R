library(igraph)
                                        # for interactive use
.plot <- function(g, ...) plot(g, vertex.size = 4, vertex.label = "", ...)

### deterministic
periodiclattice <- make_lattice(length = 10, dim = 2, circular = TRUE)
squarelattice <- make_lattice(length = 10, dim = 2, circular = FALSE)
tree_3 <- make_tree(100, 3, "undirected")
tree_4 <- make_tree(100, 4, "undirected")
tree_5 <- make_tree(100, 5, "undirected")

### stochastic
N <- 100
M <- 200
                                        # small world
sw_3 <- sample_smallworld(dim = 1, size = N, nei = 3, p = 0.05)
sw_4 <- sample_smallworld(dim = 1, size = N, nei = 4, p = 0.05)
sw_5 <- sample_smallworld(dim = 1, size = N, nei = 5, p = 0.05)
                                        # gnm
er <- sample_gnm(N, M)
                                        # geometric
ggt <- sample_grg(N, .2, TRUE)
ggs <- sample_grg(N, .2, FALSE)
                                        # ba
ba <- sample_pa(N, m= 2, directed = FALSE, start.graph = make_full_graph(3))
                                        # forestfire
ff <- sample_forestfire(N, .37, .32/.37, directed = FALSE)

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
## canton  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Cantontxt.txt", sep = "\t"))
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
## powder  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Powdertxt.txt", sep = "\t"))
stony  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Stonytxt.txt", sep = "\t"))
suttonAu  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/SuttonAutxt.txt", sep = "\t"))
suttonSp  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/SuttonSptxt.txt", sep = "\t"))
suttonSu  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/SuttonSutxt.txt", sep = "\t"))
troy  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Troytxt.txt", sep = "\t"))
venlaw  <- graph_from_biadjacency_matrix(read.table("../networks/originals/thomps_towns_text/Venlawtxt.txt", sep = "\t"))

## make graphlist
graphlist <- list(
    pl = periodiclattice, sl = squarelattice, t3 = tree_3, t4 = tree_4, t5 = tree_5,
    sw3 = sw_3, sw4 = sw_4, sw5 = sw_5, er = er, ggt = ggt, ggs = ggs, ba = ba, ff = ff,
    dolphin = dolphin, proximity = proximity, genefusion = genefusion, jazz = jazz, netsci = netsci,
    pdzbase = pdzbase, chesapeake = chesapeake,
    akatoreA = akatoreA, akatoreB = akatoreB, berwick = berwick, blackrock = blackrock, catlins = catlins,
    coweeta17 = coweeta17, coweeta1 = coweeta1, dempstersA = dempstersA, dempstersS = dempstersS,
    dempstersSu = dempstersSu, german = german, healy = healy, kyeburn = kyeburn, lilkyeburn = lilkyeburn,
    martins = martins, narrowdale = narrowdale, northcol = northcol, stony = stony, suttonAu = suttonAu,
    suttonSp = suttonSp, suttonSu = suttonSu, troy = troy, venlaw = venlaw
)
## then loop through, replacing each with simplify(largest_component(g))
for(i in seq_along(graphlist)) graphlist[[i]] <- simplify(largest_component(graphlist[[i]]))

## for checking
i <- 1

.plot(graphlist[[i]], main = names(graphlist)[i])
i <- i + 1

saveRDS(graphlist, file = "networks.rds")

## ## dolphin
## ## celegans
## celegans <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/celegans.csv")), FALSE)))
## ## proximity
## ## euroroad
## euroroad <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/euro-road.csv")), FALSE)))
## ## email
## email <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/arenas-email.csv")), FALSE)))
## ## FAA?
## faa <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/FAA.csv")), FALSE)))
## ## gene fusion?
## ## hamsterster?
## hamsterster <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/hamsterster.csv")), FALSE)))
## ## jazz?
## ## netsci?
## ## pdz base?
## ## powergrid?
## powergrid <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/powergrid.csv")), FALSE)))
## ## yeast metabolic?
## yeast <- simplify(largest_component(graph_from_edgelist(as.matrix(read.csv("../networks/yeast-metabolic.csv")), FALSE)))

## gl <- list(
##     dolphin = dolphin, celegans = celegans, proximity = proximity, euroroad = euroroad, email = email, faa = faa,
##     genefusion = genefusion, hamsterster = hamsterster, jazz = jazz, netsci = netsci, pdzbase = pdzbase,
##     ## powergrid = powergrid,
##     yeast = yeast
## )

## sapply(gl, is_simple)
## sapply(gl, is_connected)
## sapply(gl, is_directed)
## sapply(gl, vcount) # remove powergrid as too big?
## sapply(gl, ecount) # the email network has the most edges besides powergrid

## N <- as.integer(round(median(sapply(gl, vcount))))
## M <- as.integer(round(median(sapply(gl, ecount))))
          

## ### Make random graphs with N = mean(vcount(graphs above here)).
## ### definitely use ER and BA. Can also use small world and stochastic block model. Or a geometric random graph. 
## ## er
## er <- simplify(largest_component(sample_gnm(N, M)))
## ## ba
## ba <- simplify(largest_component(sample_pa(N, m = 2, directed = FALSE, start.graph = make_full_graph(3))))
## ## sw
## sw <- simplify(largest_component(sample_smallworld(dim = 1, size = N, nei = 5, p = 0.05)))
## ## gg
## gg <- simplify(largest_component(sample_grg(N, 0.1, TRUE)))

## ## sbm
## prefmat <- matrix(0.005, ncol = 6, nrow = 6)
## diag(prefmat) <- 0.05
## sizes <- rbeta(6, 2, 2)
## sizes <- round(N*(sizes/sum(sizes)))
## sbm <- simplify(largest_component(sample_sbm(N, prefmat, sizes)))
## cl <- cluster_louvain(sbm)
## V(sbm)$color <- cl$membership



## rgl <- list(er = er, ba = ba, sw = sw, gg = gg, sbm = sbm)

## graphlist <- c(gl, rgl)


