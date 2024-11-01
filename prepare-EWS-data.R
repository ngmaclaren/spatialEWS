library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sfsmisc)
library(moments)
library(sdn)
source("calc-functions.R")
palette("Okabe-Ito")
save_plots <- FALSE # TRUE

                                        # How many samples to use for classification?
n <- 5

networks <- readRDS("./data/networks.rds")

simdir <- "./data/sims/"
simfiles <- list.files(simdir, pattern = ".rds")
simdata <- lapply(paste0(simdir, simfiles), function(simfile) readRDS(simfile))

dyns <- sapply(simdata, attr, "model")
nets <- sapply(simdata, attr, "network")
bparams <- sapply(simdata, attr, "bparam")
directions <- sapply(simdata, attr, "direction")
cparam.vals <- lapply(simdata, attr, "bparam.vals")
As <- lapply(simdata, function(dat) attr(dat, "params")$A)

dfs <- lapply(simdata, promote_df)

dfnames <- mapply(
    function(dyn, net, bp, dir) paste(dyn, net, bp, dir, sep = "-"),
    dyns, nets, bparams, directions,
    USE.NAMES = FALSE
)
names(dfs) <- dfnames

conds <- as.data.frame(cbind(dyns, bparams, directions, nets))
plotorder <- with(conds, order(dyns, bparams, directions, nets))

moranIs<- mcmapply(
    function(df, A) apply(df, 1, global_moran, A = A),
    dfs, As, SIMPLIFY = FALSE, USE.NAMES = TRUE, mc.cores = ncores
)
ssds <- lapply(dfs, function(df) apply(df, 1, sd))
vars <- lapply(dfs, function(df) apply(df, 1, moment, order = 2, central = TRUE))
cvs <- lapply(dfs, function(df) apply(df, 1, CV))
skews <- lapply(dfs, function(df) {
    g_1 <- apply(df, 1, skewness)
                                        # Accounts for difference in expected change in skewness based on the
                                        # direction of the simulations:
    if(attr(df, "direction") == "down") g_1 <- -g_1
    g_1
})
kurts <- lapply(dfs, function(df) apply(df, 1, kurtosis))
taus <- list(
    moranI = mapply(get_tau, dfs, moranIs, TRUE),
    ssd = mapply(get_tau, dfs, ssds, TRUE),
    var = mapply(get_tau, dfs, vars, TRUE),
    cv = mapply(get_tau, dfs, cvs, TRUE),
    skew = mapply(get_tau, dfs, skews, TRUE),
    kurt = mapply(get_tau, dfs, kurts, TRUE)
)

save.image("./data/EWS-data.RData")
