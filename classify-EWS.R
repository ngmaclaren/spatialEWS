library(igraph)
library(localsolver)
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
plotorder <- with(conds, order(dyns, bparams, directions, nets))

u.down <- which(grepl("-u-", dfnames, fixed = TRUE) & grepl("-down", dfnames, fixed = TRUE))
u.up <- which(grepl("-u-", dfnames, fixed = TRUE) & grepl("-up", dfnames, fixed = TRUE))
D.down <- which(grepl("-D-", dfnames, fixed = TRUE) & grepl("-down", dfnames, fixed = TRUE))
D.up <- which(grepl("-D-", dfnames, fixed = TRUE) & grepl("-up", dfnames, fixed = TRUE))

moranIs<- mapply(
    function(df, A) apply(df, 1, global_moran, A = A),
    dfs, As, SIMPLIFY = FALSE, USE.NAMES = TRUE
)
ssds <- lapply(dfs, function(df) apply(df, 1, sd))

basins <- list( ## GLOBAL
    doublewell = 3,
    genereg = 0.1,
    mutualistic = 1,
    SIS = 0.1
)

## Now, all the information I need should be with the df
get_idx <- function(df) {
    direction <- attr(df, "direction")
    dynamics <- attr(df, "model")
    basin <- basins[[dynamics]]

    switch(
        direction,
        down = which(apply(df, 1, min) > basin),
        up = which(apply(df, 1, max) < basin)
    )
}

get_samples <- function(df, which = c("near", "far"), n = 3) {
    whichslope <- match.arg(which)

    idx <- get_idx(df)

    switch(
        whichslope,
        near = rev(seq(length(idx), by = -1, length.out = n)),
        far = seq(1, by = 1, length.out = n)
    )
}

get_slope <- function(df, EWS, which = c("near", "far"), n = 3, return.model = FALSE) {
    whichslope <- match.arg(which)

    samples <- get_samples(df, whichslope, n)

    cparam <- attr(df, "bparam.vals")[samples]
    ews <- EWS[samples]
    m <- lm(ews ~ cparam)

    if(return.model) {
        return(m)
    } else {
        return(as.numeric(coef(m)[2]))
    }
}

                                        # test
slopes <- data.frame(
    ssd.far = mapply(get_slope, dfs, ssds, which = "far", n = 3),
    ssd.near = mapply(get_slope, dfs, ssds, which = "near", n = 3),
    moran.far = mapply(get_slope, dfs, moranIs, which = "far", n = 3),
    moran.near = mapply(get_slope, dfs, moranIs, which = "near", n = 3)
)

                                        # I think we don't want an absolute value threshold: no way to know in
                                        # advance what scale the slope differences will be on.
classify <- function(df, EWS, multiple = 2) {#, threshold = 2) {
    firstslope <- get_slope(df, EWS, "far")
    secondslope <- get_slope(df, EWS, "near")

    firstsign <- sign(firstslope)
    secondsign <- sign(secondslope)
    mult <- secondslope/firstslope

    direction <- attr(df, "direction")

    if(direction == "up") {
        if(sign(firstslope) > 0 & sign(secondslope) > 0 & abs(mult) > multiple) {# & secondslope > threshold) {
            return("consistent")
        } else if(sign(firstslope) < 0 & sign(secondslope) > 0 & abs(mult) > multiple) {# & secondslope > threshold) {
            return("signswitch")
        } else {
            return("bad")
        }
    } else { # down
        if(sign(firstslope) < 0 & sign(secondslope) < 0 & abs(mult) > multiple) {# & secondslope > threshold) {
            return("consistent")
        } else if(sign(firstslope) > 0 & sign(secondslope) < 0 & abs(mult) > 0.5*multiple) {# & secondslope > threshold) {
            return("signswitch")
        } else {
            return("bad")
        }
    }
}

diagnostic_plot <- function(df, moranI, ssd, main = "") {
    cparam.vals <- attr(df, "bparam.vals")
    idx <- get_idx(df)
    farsample <- get_samples(df, "far")
    nearsample <- get_samples(df, "near")

    farmodel.moran <- get_slope(df, moranI, "far", return.model = TRUE)
    farmodel.ssd <- get_slope(df, ssd, "far", return.model = TRUE)
    nearmodel.moran <- get_slope(df, moranI, "near", return.model = TRUE)
    nearmodel.ssd <- get_slope(df, ssd, "near", return.model = TRUE)

    test.moran <- classify(df, moranI)
    test.ssd <- classify(df, ssd)

    bifplot(df, cparam.vals, col = adjustcolor(1, 0.5), lwd = 0.5, main = main)
    mtext("Test results", line = 2, adj = 0, font = 2)
    mtext(paste("Moran's I:", test.moran), line = 1, adj = 0)
    mtext(paste("Spatial SD:", test.ssd), line = 0, adj = 0)
    mtext(paste("Dynamics:", attr(df, "model")), line = 2, adj = 1)
    mtext(paste("Control parameter:", attr(df, "bparam")), line = 1, adj = 1)
    mtext(paste("Direction:", attr(df, "direction")), line = 0, adj = 1)

    par(new = TRUE)
    plot(
        cparam.vals[idx], moranI[idx], type = "l", col = adjustcolor(2, 0.5), lwd = 2, lty = 3,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    lines(cparam.vals[farsample], predict(farmodel.moran), lwd = 3, lty = 1, col = 2)
    lines(cparam.vals[nearsample], predict(nearmodel.moran), lwd = 3, lty = 1, col = 2)

    par(new = TRUE)
    plot(
        cparam.vals[idx], ssd[idx], type = "l", col = adjustcolor(3, 0.5), lwd = 2, lty = 3,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    lines(cparam.vals[farsample], predict(farmodel.ssd), lwd = 3, lty = 1, col = 3)
    lines(cparam.vals[nearsample], predict(nearmodel.ssd), lwd = 3, lty = 1, col = 3)

    legend("bottomright", bty = "n", lwd = 2, col = 2:3, lty = 1,
           legend = c("Moran's I", "Spatial std. dev."))
}

pdf("./img/diagnostic-plots.pdf")
mapply(diagnostic_plot, dfs[plotorder], moranIs[plotorder], ssds[plotorder])
dev.off()

results <- data.frame(
    result = c(mapply(classify, dfs, moranIs), mapply(classify, dfs, ssds)),
    EWS = c(rep("moranI", length(dfs)), rep("ssd", length(dfs))),
    dynamics = dyns,
    network = nets,
    bparam = bparams,
    direction = directions
)

restab <- table(results$EWS, results$result)
proportions(restab)*2

table(results$dynamics, results$result, results$EWS)

table(results$direction, results$result, results$EWS)
