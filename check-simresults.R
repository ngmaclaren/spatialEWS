                                        # File for examining the results of the simulations. For debugging,
                                        # determining thresholds, etc.

library(igraph)
library(sdn)

networks <- readRDS("./data/networks.rds")
netlist <- scan("./data/newnetworknames.txt", character())
## networks <- networks[nets]

simdir <- "./data/sims/"
simfiles <- list.files(simdir, pattern = ".rds")
simdata <- lapply(paste0(simdir, simfiles), function(simfile) readRDS(simfile))

dyns <- sapply(simdata, attr, "model")
nets <- sapply(simdata, attr, "network")
bparams <- sapply(simdata, attr, "bparam")
directions <- sapply(simdata, attr, "direction")
cparam.vals <- lapply(simdata, attr, "bparam.vals")
As <- lapply(simdata, function(dat) attr(dat, "params")$A)

dfs <- lapply(
    simdata,
    function(sim) {
        attribs <- attributes(sim)
        df <- sim[[1]]
        attributes(df) <- c(attributes(df), attribs)
        df
    }
)

dfnames <- mapply(
    function(dyn, net, bp, dir) paste(dyn, net, bp, dir, sep = "-"),
    dyns, nets, bparams, directions,
    USE.NAMES = FALSE
)
names(dfs) <- dfnames

conds <- as.data.frame(cbind(dyns, bparams, directions, nets))
plotorder <- with(conds, order(dyns, bparams, directions, nets))#[grep("genereg", dfnames)]

pdf("./img/checksims.pdf")
mapply(
    function(df, cparam.val, main) {
        bifplot(df, cparam.val, lwd = 0.5, col = 1, main = main)#, ylim = c(0, 0.1))
        ##abline(h = 0.001, col = "red")
        ##print(main)
    }, 
    dfs[plotorder], cparam.vals[plotorder], dfnames[plotorder]
)
dev.off()
