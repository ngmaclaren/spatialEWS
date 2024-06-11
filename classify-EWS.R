library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
source("calc-functions.R")
palette("Okabe-Ito")
## with(list(n = length(palette())), plot(seq(n), rep(1, n), col = seq(n), cex = 5, pch = 16))
## dev.new(width = 14, height = 3); plot(0:24, rep(1, 25), pch = 0:24, cex = 5, col = 1, bg = 2)
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
##     function(sim) {
##         attribs <- attributes(sim)
##         df <- sim[[1]]
##         attributes(df) <- c(attributes(df), attribs)
##         df
##     }
## )

dfnames <- mapply(
    function(dyn, net, bp, dir) paste(dyn, net, bp, dir, sep = "-"),
    dyns, nets, bparams, directions,
    USE.NAMES = FALSE
)
names(dfs) <- dfnames

conds <- as.data.frame(cbind(dyns, bparams, directions, nets))
plotorder <- with(conds, order(dyns, bparams, directions, nets)) # makes it alphabetical. Should specify. Want the double-well type and the transcritical type showing up together.

## u.down <- which(grepl("-u-", dfnames, fixed = TRUE) & grepl("-down", dfnames, fixed = TRUE))
## u.up <- which(grepl("-u-", dfnames, fixed = TRUE) & grepl("-up", dfnames, fixed = TRUE))
## D.down <- which(grepl("-D-", dfnames, fixed = TRUE) & grepl("-down", dfnames, fixed = TRUE))
## D.up <- which(grepl("-D-", dfnames, fixed = TRUE) & grepl("-up", dfnames, fixed = TRUE))

### Spectral reddening as EWS? ###
## df <- dfs[[1]]
## dl <- apply(df, 1, stats::spectrum, plot = FALSE)
## test <- sapply(dl, `[[`, "spec")
## rownames(test) <- dl[[1]]$freq
## colnames(test) <- cparam.vals[[1]]
## freqs <- c(0, dl[[1]]$freq)
## Ds <- cparam.vals[[1]]
## Ds <- c(Ds, 0.3-0.007)
## image(x = rev(Ds), y = freqs, z = t(test), xlab = "D", ylab = "Frequency")

moranIs<- mcmapply(
    function(df, A) apply(df, 1, global_moran, A = A),
    dfs, As, SIMPLIFY = FALSE, USE.NAMES = TRUE, mc.cores = ncores
)
ssds <- lapply(dfs, function(df) apply(df, 1, sd))

                                        # Restart here if adjusting basins                                #####
taus <- list(
    moranI = mapply(get_tau, dfs, moranIs, TRUE),
    ssd = mapply(get_tau, dfs, ssds, TRUE)
)
## slopes <- data.frame(
##     ssd.far = mapply(get_slope, dfs, ssds, which = "far", n = n),
##     ssd.near = mapply(get_slope, dfs, ssds, which = "near", n = n),
##     moran.far = mapply(get_slope, dfs, moranIs, which = "far", n = n),
##     moran.near = mapply(get_slope, dfs, moranIs, which = "near", n = n)
## )

plotdata <- data.frame(
    tau.I = taus$moranI, # these are 'sign-adjusted'
    tau.sd = taus$ssd,
    ## slope.far.I = slopes$moran.far,
    ## slope.near.I = slopes$moran.near,
    ## slope.far.sd = slopes$ssd.far,
    ## slope.near.sd = slopes$ssd.near,
    dynamics = dyns,
    network = nets,
    cparam = bparams,
    direction = directions,
    row.names = seq_along(dfs)
)
plotdata$dynamics <- factor(plotdata$dynamics, levels = c("doublewell", "mutualistic", "SIS", "genereg"))
plotdata$network <- factor(plotdata$network)
plotdata$cparam <- factor(plotdata$cparam)
plotdata$direction <- factor(plotdata$direction)

## then maybe a strip plot by qualitative var
## Color by EWS  #### START HERE: separate more. Within a dynamics, four rows: I D, SD D, I u, SD u
ht <- 5; wd <- 12
if(save_plots) {
    pdf("./img/tauplot.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
labelsize <- 1.75
ptcex <- 2
par(mfrow = c(1, 2), mar = c(4, 10, 0.5, 0.5))
with(subset(plotdata, direction == "down"), {
    xlim <- c(-1, 1)
    ylim <- rev(c(0.5, 4.5))
    pchs <- with(list(x = as.numeric(cparam)), ifelse(x == 1, x, 4))
    ##adj <- as.numeric(cparam)/4
    plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    abline(v = 0)
    axis(2, at = 1:4, labels = levels(dynamics), cex.axis = labelsize, las = 2)
    axis(1, cex.axis = labelsize)
    title(xlab = "Sign-adjusted tau", cex.lab = labelsize)
    mtext("(a)", line = -1.5, adj = 0.01, cex = labelsize)
    box()
    points(tau.I, jitter(as.numeric(dynamics) - 0.2, amount = 0.1), pch = pchs, col = 2, cex = ptcex)
    points(tau.sd, jitter(as.numeric(dynamics) + 0.2, amount = 0.1), pch = pchs, col = 3, cex = ptcex)
})
with(subset(plotdata, direction == "up"), {
    xlim <- c(-1, 1)
    ylim <- rev(c(0.5, 4.5))
    pchs <- with(list(x = as.numeric(cparam)), ifelse(x == 1, x, 4))
    plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    abline(v = 0)
    axis(2, at = 1:4, labels = levels(dynamics), cex.axis = labelsize, las = 2)
    axis(1, cex.axis = labelsize)
    title(xlab = "Sign-adjusted tau", cex.lab = labelsize)
    mtext("(b)", line = -1.5, adj = 0.01, cex = labelsize)
    box()
    points(tau.I, jitter(as.numeric(dynamics) - 0.2, amount = 0.1), pch = pchs, col = 2, cex = ptcex)
    points(tau.sd, jitter(as.numeric(dynamics) + 0.2, amount = 0.1), pch = pchs, col = 3, cex = ptcex)
    legend("bottomright", bty = "n", ncol = 1, col = c(2:3, 9, 9), pch = c(1, 1, 1, 4),
           cex = 0.75*labelsize, pt.cex = 0.75*ptcex,
           legend = c("Moran's I", "SSD", "D", "u"))
})
if(save_plots) dev.off()

                                        # I think we don't want an absolute value threshold: no way to know in
                                        # advance what scale the slope differences will be on.

diagnostic_plot <- function(df, moranI, ssd, n = 3, main = "", ...) {
    cparam.vals <- attr(df, "bparam.vals")
    idx <- get_idx(df)
    farsample <- get_samples(df, "far", n = n)
    nearsample <- get_samples(df, "near", n = n)

    farmodel.moran <- get_slope(df, moranI, "far", n = n, return.model = TRUE)
    farmodel.ssd <- get_slope(df, ssd, "far", n = n, return.model = TRUE)
    nearmodel.moran <- get_slope(df, moranI, "near", n = n, return.model = TRUE)
    nearmodel.ssd <- get_slope(df, ssd, "near", n = n, return.model = TRUE)

    test.moran <- classify(df, moranI, n = n)
    test.ssd <- classify(df, ssd, n = n)

    bifplot(df, cparam.vals, col = adjustcolor(1, 0.5), lwd = 0.5, ylim = c(0, 0.02), ...)

    mtext("Test results", line = 2, adj = 0, font = 2)
    mtext(paste("Moran's I:", test.moran), line = 1, adj = 0)
    mtext(paste("Spatial SD:", test.ssd), line = 0, adj = 0)

    mtext("Kendall's tau", line = 2, adj = 0.5, font = 2)
    mtext(paste("Moran's I:", round(get_tau(df, moranI, TRUE), 2)), line = 1, adj = 0.5)
    mtext(paste("Spatial SD:", round(get_tau(df, ssd, TRUE), 2)), line = 0, adj = 0.5)

    mtext(paste("Dynamics:", attr(df, "model")), line = 3, adj = 1)
    mtext(paste("Control parameter:", attr(df, "bparam")), line = 2, adj = 1)
    mtext(paste("Direction:", attr(df, "direction")), line = 1, adj = 1)
    mtext(paste("Network:", attr(df, "network")), line = 0, adj = 1)

    par(new = TRUE)
    plot(
        cparam.vals[idx], moranI[idx], type = "l", col = 2, lwd = 3, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    lines(cparam.vals[farsample], predict(farmodel.moran), lwd = 1.5, lty = 1, col = 1)
    lines(cparam.vals[nearsample], predict(nearmodel.moran), lwd = 1.5, lty = 1, col = 1)

    par(new = TRUE)
    plot(
        cparam.vals[idx], ssd[idx], type = "l", col = 3, lwd = 3, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    lines(cparam.vals[farsample], predict(farmodel.ssd), lwd = 1.5, lty = 1, col = 1)
    lines(cparam.vals[nearsample], predict(nearmodel.ssd), lwd = 1.5, lty = 1, col = 1)

    legend("bottomright", bty = "n", lwd = 2, col = 2:3, lty = 1,
           legend = c("Moran's I", "Spatial std. dev."))
}

pdf("./img/diagnostic-plots.pdf")
## mapply(diagnostic_plot, dfs[plotorder], moranIs[plotorder], ssds[plotorder], n = n)
with(
    list(theseplots = grep("genereg", dfnames)), 
    mapply(diagnostic_plot, dfs[theseplots], moranIs[theseplots], ssds[theseplots], n = n)
)
dev.off()

results <- data.frame(
    result = c(mapply(classify, dfs, moranIs), mapply(classify, dfs, ssds)),
    EWS = c(rep("moranI", length(dfs)), rep("ssd", length(dfs))),
    dynamics = dyns,
    network = nets,
    bparam = bparams,
    direction = directions
)

## restab <- table(results$EWS, results$result)
## proportions(restab)*2
## addmargins(table(results$dynamics, results$result, results$EWS))[, , 1:2]
## addmargins(table(results$direction, results$result, results$EWS))[, , 1:2]
## table(results$dynamics, results$direction,  results$result, results$EWS)

print("Moran's I")
with(list(df = subset(results, EWS == "moranI")), table(df$dynamics, df$bparam, df$result, df$direction))
print("Spatial standard deviation")
with(list(df = subset(results, EWS == "ssd")), table(df$dynamics, df$bparam, df$result, df$direction))


                                        # data for example figure panels
with(
    list(
        idx = grep("drug", nets)
    ), {
        print(taus$ssd[idx])
        print(mapply(classify, dfs[idx], ssds[idx]))
    }
)
