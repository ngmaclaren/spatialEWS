## Example network is the drug interaction network.
## Use layout with specific (i.e., in in/cm) margins, etc., to make all plots exactly the same size.
## Six panels:
## A: dw, u, up
## B: dw, u, down
## C: dw, D, up
## D: dw, D, down
## E: SIS, D, up
## G: SIS, D, down
## Otherwise, each plot should be basically the same. The plotit() function should work? Let's try it.

library(sfsmisc)
library(sdn)
source("calc-functions.R")
save_plots <- TRUE # FALSE
palette("Okabe-Ito")

plotit <- function(sim, col.xi = 9, col.ssd = 3, pad.ssd = NULL, labelsize = NULL, ...) {
                                        # assumes palette is Okabe-Ito
    X <- promote_df(sim)
    ##X <- sim[[1]]
    bparam.vals <- attr(X, "bparam.vals")
    direction <- attr(X, "direction")
    dynamics <- attr(X, "model")
    basin.limit <- get(paste0(".", dynamics))$basin.limit # b/c fixed since ran sims

                                        # Pick the bparam vals at which to compute ssd
    edges <- switch(
        direction,
        down = apply(X, 1, min),
        up = apply(X, 1, max)
    )
    far1 <- 10
    far2 <- 30
    nearest <- switch(
        direction,
        down = which(edges < basin.limit)[1] - 1,
        up = which(edges > basin.limit)[1] - 1
    )
    near <- nearest - 3
        
                                        # data for marking lines, as matrix
    markvec <- c(far1, far2, near, nearest)
    mark.bpv <- bparam.vals[markvec]
    mark.X <- X[markvec, ]
    if(is.null(pad.ssd)) pad.ssd <- .01*mean(mark.X)
                                        # !! Place Moran's I code here if showing

                                        # Plot
    if(is.null(labelsize)) labelsize <- 2
    par(mar = c(5, 5, 0.5, 0.5))
    matplot(
        bparam.vals, X, type = "l", lty = 1, lwd = 0.5, col = col.xi,
        xlab = "", ylab = "", axes = FALSE, cex.lab = 2,
        ...
    )
    ##lines(bparam.vals, rowMeans(X), lty = 1, lwd = 4, col = 5)
    eaxis(1, cex.axis = 2)
    eaxis(2, cex.axis = 2)
    title(xlab = attr(sim, "bparam"), cex.lab = labelsize)
    title(ylab = expression(x[i]), cex.lab = labelsize, line = 3.2)
    segments(
        x0 = mark.bpv, y0 = apply(mark.X, 1, min) - pad.ssd, y1 = apply(mark.X, 1, max) + pad.ssd,
        col = col.ssd, lty = 1, lwd = 4
    )

    sds <- apply(mark.X, 1, sd)
    cvs <- apply(mark.X, 1, CV)
    skews <- apply(mark.X, 1, moments::skewness)
    if(attr(X, "direction") == "down") skews <- -skews

    allsd <- apply(X, 1, sd)
    allcv <- apply(X, 1, CV)
    allskew <- apply(X, 1, moments::skewness)
    if(attr(X, "direction") == "down") allskew <- -allskew
    ## mtext(tools::toTitleCase(classify(X, allsd, n = 5)), col = col.ssd, line = -1)
    mtext(tools::toTitleCase(classify(X, allcv, n = 5)), col = col.ssd, line = -1)
    mtext(tools::toTitleCase(classify(X, allskew, n = 5)), col = col.ssd+1, line = -3)

    
    cat("\n", paste(dynamics, direction, attr(X, "bparam")),
        "\nControl parameter values:\n", round(mark.bpv, 3),
        ##"\nSpatial standard deviation:\n", round(sds, 3),
        "\nCoefficient of variation:\n", round(cvs, 3),
        "\nSpatial skewness:\n", round(skews, 3),
        "\n"
        )
}

sims <- list(
    dw.uu = readRDS("./data/sims/doublewell-drug-u-sde-50-up-0-1-TRUE.rds"), ## A
    dw.ud = readRDS("./data/sims/doublewell-drug-u-sde-50-down-0-1-TRUE.rds"), ## B
    dw.Du = readRDS("./data/sims/doublewell-drug-D-sde-50-up-0-1-TRUE.rds"), ## C
    dw.Dd = readRDS("./data/sims/doublewell-drug-D-sde-50-down--5-1-TRUE.rds"), ## D
    SIS.Du = readRDS("./data/sims/SIS-drug-D-sde-50-up-0-1-TRUE.rds"), ## E
    SIS.Dd = readRDS("./data/sims/SIS-drug-D-sde-50-down-0-1-TRUE.rds")## F
)
## lat.Dd <- readRDS("./data/sims/doublewell-lattice-D-sde-50-down--5-1-TRUE.rds")

pht <- pwd <- 5
wd <- 3*pwd
ht <- 2*pwd
ylims <- list(
    dw.uu = c(0.9, 2.5),
    dw.ud = c(4, 5.5),
    dw.Du = c(0.9, 2.5),
    dw.Dd = c(3.5, 8.5),
    SIS.Du = c(0, 0.05),
    SIS.Dd = c(0, 1)
)
xlims <- lapply(sims[1:4], function(sim) range(attr(sim, "bparam.vals")))
xlims$SIS.Du <- c(0, 0.15)
xlims$SIS.Dd <- c(0, 1)
    

if(save_plots) {
    pdf("./img/example-drug-input.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(2, 3))
## lapply(sims, plotit)
mapply(function(sim, ylim, xlim) plotit(sim, ylim = ylim, xlim = xlim), sims, ylims, xlims)
if(save_plots) dev.off()


moran_components <- function(x, A) {
    N <- length(x)
    stopifnot(N == nrow(A))
    stopifnot(nrow(A) == ncol(A))

    x <- x - mean(x)

    num <- sum(A*outer(x, x))
    denom <- sum(x^2)

    c(numerator = num, denominator = denom)
}

get_moran_components <- function(sim) {
    df <- promote_df(sim)
    idx <- get_idx(df)

    g <- readRDS("./data/networks.rds")$drug
    A <- igraph::as_adj(g, "both", sparse = FALSE)

    t(apply(df, 1, moran_components, A = A)[, idx])
}


## with(list(df = promote_df(sims$dw.uu)), apply(df, 1, moran_components, A = A))
## get_moran_components(sims$dw.uu)

plot_moran_components <- function(sim, labelsize = 2) {
    g <- readRDS("./data/networks.rds")$drug
    A <- igraph::as_adj(g, "both", sparse = FALSE)

    df <- promote_df(sim)

    mc <- get_moran_components(sim)
    mI <- apply(df, 1, global_moran, A = A)
    bpv <- attr(df, "bparam.vals")
    
    idx <- get_idx(df)

    Icol <- 2
    numcol <- 6
    denomcol <- 7

    par(mar = rep(4.75, 4))
    matplot(bpv[idx], mc, col = c(numcol, denomcol),
            type = "l", lty = 1, lwd = 2,
            xlim = range(bpv),
            xlab = attr(df, "bparam"), ylab = "Variance", axes = FALSE, cex.lab = labelsize)
    eaxis(1, cex.axis = 0.75*labelsize)
    eaxis(2, cex.axis = 0.75*labelsize)

    par(new = TRUE)
    plot(
        bpv[idx], mI[idx], type = "l", col = Icol, lwd = 2,
        xlim = range(bpv),# ylim = c(-0.75, 0.75),
        xlab = "", ylab = "", axes = FALSE
    )
    eaxis(4, cex.axis = 0.75*labelsize, col = Icol, col.axis = Icol, small.args = list(col = Icol))
    
    legend("topleft", bty = "n", col = c(Icol, numcol, denomcol), lty = 1, lwd = 2,
           legend = c("Moran's I", "Numerator", "Denominator"))
    mtext(tools::toTitleCase(classify(df, mI, n = 5)))
}

if(save_plots) {
    pdf("./img/example-moran-input.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(2, 3))
lapply(sims, function(sim) {
    plot_moran_components(sim)
})
if(save_plots) dev.off()
