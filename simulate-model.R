## In this version, set better parameters using a table stored separately

library(optparse)
option_list <- list(
    make_option(
        c("-m", "--model"), type = "character", default = "doublewell",
        help = "Default is %default. Available: doublewell, mutualistic, genereg, and SIS (from package:localsolver). Each has some standard parameter values and accepts /D/, a node susceptibility parameter, and /u/, a stress parameter. As a control parameter, /u/ can be positive with '--direction=up' or negative with '--direction=down'."
    ),
    make_option(
        c("-g", "--network"), type = "character", default = "dolphin",
        help = "Default is %default. Available: montreal, dolphin, erdosrenyi, smallworld, barabasialbert, tree, lattice, adjnoun, and students."
    ),
    make_option(
        c("-p", "--bparam"), type = "character", default = "D",
        help = "Default is %default. Available: D or u."
    ),
    make_option(
        c("-k", "--simkind"), type = "character", default = "sde",
        help = "Default is %default. Available: ode or sde."
    ),
    make_option(
        c("-t", "--simtime"), type = "integer", default = 50,
        help = "Default is %default. Length of simulation in user time units (not Δt)."
    ),
    make_option(
        c("-d", "--direction"), type = "character", default = "up",
        help = "For models that can go up (more positive) or down (more negative). Default is %default. Make sure to choose an appropriate initial value for x."
    ),
    make_option(
        c("-s", "--sigma"), type = "numeric", default = NaN,
        help = "The strength of the noise process, if using. Default is %default. If finite, will use in place of the standard value."
    ),
    make_option(
        c("-x", "--xinit"), type = "numeric", default = NaN,
        help = "Default is %default. If finite, will use in place of the standard value, which is set according to the direction (up or down)."
    ),
    make_option(
        ## Calling this uinit because I still want to allow an initial u value, perhaps heterogeneous, that still works with using u as the bifurcation parameter. 
        c("-u", "--uinit"), type = "numeric", default = NaN,
        help = "Default is %default. If finite, will use in place of the standard value for node stress value. The standard value is zero for all models."
    ),
    make_option(
        c("-D", "--Dinit"), type = "numeric", default = NaN,
        help = "Default is %default. If finite, will use in place of the standard value for node susceptibility (which we also refer to as connection strength: when D is a scalar, the two concepts can be written the same way). The standard value depends on the model."
    ),
    make_option(
        c("-n", "--ntrials"), type = "integer", default = 1,
        help = "The number of independent simulations to run; default is %default. Output is a list of data frames of length `ntrials`, the values of which are equilibrium node states for each node (columns) and bifurcation parameter (rows). There will be additional information in the object attributes."
    )
)

args <- parse_args(
    OptionParser(option_list = option_list),
    convert_hyphens_to_underscores = TRUE
)

                                        # For debugging
if(interactive()) {
    args$network <- "adjnoun"
    args$model <- "doublewell"
    args$bparam <- "u"
    args$direction <- "down"
    ## args$uinit <- -5
}

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(localsolver)

networks <- readRDS("./data/networks.rds")

## order is: dyn, net, bparam, simkind, simtime, direction, ntrials
## savedir <- "/projects/academic/naokimas/neil/spatialEWS/"
savedir <- "./data/sims/"
filetag <- paste(unlist(args)[-which(unlist(args) %in% c("NaN", "FALSE"))], collapse = "-")
print(filetag)

## source("functions.R")
## load(paste0("./data/", args$network, ".rda")) # should convert to .rds at some point

## g <- get(args$network)
g <- networks[[args$network]]
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)

deltaT <- 0.01
lmax <- 100

model <- get(args$model)
modelparams <- get(paste0(".", args$model))

                                        # Set params from the separate file here
refparams <- read.csv("./data/cparam-values.csv")
rowfilter <- refparams$model == args$model & refparams$cparam == args$bparam & refparams$direction == args$direction
farlim <- refparams$farlim[rowfilter]
biflim <- refparams$biflim[rowfilter]
if(args$bparam == "D") {
    modelparams$Ds <- seq(farlim, biflim, length.out = lmax)
} else {
    modelparams$us <- seq(farlim, biflim, length.out = lmax)
}

                                        # The initial x value often has a value that makes sense.
                                        # That value is stored in the standard parameter list.
                                        # Select that value, then assign it to every node.
if(is.finite(args$xinit)) {
    ## xinit <- rep(args$xinit, N)
    xinit <- args$xinit
} else {
    xinit <- switch(args$direction, up = modelparams$xinit.low, down = modelparams$xinit.high)
}
xinit <- rep(xinit, N)
                                        # For parameters, modify directly before including A
if(is.finite(args$sigma)) modelparams$sigma <- args$sigma
if(is.finite(args$uinit)) modelparams$u <- args$uinit
if(is.finite(args$Dinit)) modelparams$D <- args$Dinit

params <- c(modelparams, list(A = A))
control <- list(ncores = ncores, times = 0:args$simtime, deltaT = deltaT)

rng <- switch(
    args$bparam,
    ## D = switch(args$direction, up = modelparams$Ds.up, down = modelparams$Ds.down),
    D = params$Ds,
    ## u = switch(args$direction, up = modelparams$us.up, down = modelparams$us.down)
    u = params$us    
)

                                        # Debugging
if(interactive()) { # system.time() won't print to stdout like this---run each line individually, or paste somewhere else
    system.time(X <- sde(xinit, control$times, model, params, control))
    system.time(result <- solve_in_range(rng, args$bparam, model, xinit, params, control, kind = args$simkind))
    rowMeans(result)
    ## result <- solve_in_range(rng, args$bparam, model, xinit, params, control, kind = args$simkind)
    pdf(paste0(filetag, "-bifplot.pdf"))
    bifplot(result, rng, col = 1)
    abline(h = 0, col = 2)
    dev.off()
}


                                        # Main simulations
                                        # want `result` to always be a list for consistency
                                        # seeds <- ... and lapply ensure different seeds for each sim
seeds <- sample(1:100000, args$ntrials)
result <- lapply(
    seeds, function(seed) {
        set.seed(seed)
        solve_in_range(rng, args$bparam, model, xinit, params, control, kind = args$simkind)
    }
)

                                        # assign attributes
attr(result, "model") <- args$model
attr(result, "network") <- args$network
attr(result, "simkind") <- args$simkind
attr(result, "simtime") <- args$simtime
attr(result, "bparam") <- args$bparam
attr(result, "direction") <- args$direction
attr(result, "bparam.vals") <- rng
attr(result, "params") <- params
attr(result, "control") <- control
attr(result, "ntrials") <- args$ntrials

saveRDS(result, file = paste0(savedir, filetag, ".rds"))
