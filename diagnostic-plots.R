library(igraph)
library(sdn)
source("calc-functions.R")
palette("Okabe-Ito")

networks <- readRDS("./data/networks.rds")

simdir <- "./data/sims/"
simfiles <- list.files(simdir)
simdata <- lapply(paste0(simdir, simfiles), function(simfile) readRDS(simfile))
dyns <- sapply(simdata, attr, "model")
nets <- sapply(simdata, attr, "network")
bparams <- sapply(simdata, attr, "bparam")
directions <- sapply(simdata, attr, "direction")

cparam.vals <- lapply(simdata, attr, "bparam.vals")
As <- lapply(simdata, function(dat) attr(dat, "params")$A)

dfs <- lapply(simdata, `[[`, 1) # missing the down direction
dfnames <- mapply(
    function(dyn, net, bp, dir) paste(dyn, net, bp, dir, sep = "-"),
    dyns, nets, bparams, directions,
    USE.NAMES = FALSE
)
names(dfs) <- dfnames

## As <- lapply(nets, function(net) as_adj(networks[[net]], "both", sparse = FALSE))

moranIs<- mapply(
    function(df, A) apply(df, 1, global_moran, A = A),
    dfs, As, SIMPLIFY = FALSE, USE.NAMES = TRUE
)
cssds <- lapply(dfs, function(df) apply(df, 1, sd))

basins <- list(
    doublewell = 3,
    genereg = 0.1,
    mutualistic = 1,
    SIS = 0.1
)

pdf("./img/diagnostic-plots.pdf")

for(i in seq_along(dfs)) {

    bifplot(dfs[[i]], cparam.vals[[i]], col = adjustcolor(1, 0.5), lwd = 0.5, main = dfnames[i])

    basin <- basins[[dyns[i]]]
    idx <- switch(
        directions[i],
        down = which(apply(dfs[[i]], 1, min) > basin),
        up = which(apply(dfs[[i]], 1, max) < basin)
    )

    par(new = TRUE)
    plot(
        cparam.vals[[i]][idx], moranIs[[i]][idx], type = "l", col = 2, lwd = 2, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals[[i]])
    )
    par(new = TRUE)
    plot(
        cparam.vals[[i]][idx], cssds[[i]][idx], type = "l", col = 3, lwd = 2, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals[[i]])
    )

    if(i == 1) legend("bottomright", bty = "n", lwd = 2, col = 2:3, lty = 1,
                      legend = c("Moran's I", "Cross-sectional SD"))
    
}

dev.off()
