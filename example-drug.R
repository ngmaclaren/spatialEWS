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
    eaxis(1, cex.axis = 2)
    eaxis(2, cex.axis = 2)
    title(xlab = attr(sim, "bparam"), cex.lab = labelsize)
    title(ylab = expression(x[i]), cex.lab = labelsize, line = 3.2)
    segments(
        x0 = mark.bpv, y0 = apply(mark.X, 1, min) - pad.ssd, y1 = apply(mark.X, 1, max) + pad.ssd,
        col = col.ssd, lty = 1, lwd = 4
    )

    cat("\n", paste(dynamics, direction, attr(X, "bparam")),
        "\nControl parameter values:\n", mark.bpv,
        "\nSpatial standard deviation:\n", apply(mark.X, 1, sd), "\n"
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
    pdf("./img/example-drug.pdf", height = ht, width = wd)
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

with(list(df = promote_df(sims$dw.uu)), apply(df, 1, moran_components, A = A))

get_moran_components(sims$dw.uu)

plot_moran_components <- function(sim) {
    df <- promote_df(sim)
    mc <- get_moran_components(sim)
    bpv <- attr(df, "bparam.vals")
    idx <- get_idx(df)

    matplot(bpv[idx], mc, col = c("#f6d32d", "#9141ac"),
            type = "l", lty = 1, lwd = 2,
            xlim = range(bpv),
            xlab = "Control parameter", ylab = "")
    legend("bottomright", bty = "n", col = c("#f6d32d", "#9141ac"), lty = 1, lwd = 2,
           legend = c("Numerator", "Denominator"))
}

moranIs <- lapply(sims, function(sim) {
    g <- readRDS("./data/networks.rds")$drug
    A <- igraph::as_adj(g, "both", sparse = FALSE)
    df <- promote_df(sim)
    apply(df, 1, global_moran, A = A)
})
ssds <- lapply(sims, function(sim) {
    df <- promote_df(sim)
    apply(df, 1, sd)
})

mapply(function(sim, EWS) classify(promote_df(sim), EWS, n = 5), sims, moranIs)
mapply(function(sim, EWS) classify(promote_df(sim), EWS, n = 5), sims, ssds)

lapply(sims, function(sim) {
    dev.new()
    plot_moran_components(sim)
})

## pht <- pwd <- 5
## wd <- 3*pwd
## ht <- 2*pwd
## if(save_plots) {
##     pdf("./img/example-dw-v3.pdf", height = ht, width = wd)
## } else {
##     dev.new(height = ht, width = wd)
## }
## lyt <- layout(
##     matrix(c(
##         1, 1, 2, 2, 3, 3,
##         1, 1, 2, 2, 3, 3,
##         4, 4, 5, 5, 6, 6,
##         4, 4, 5, 5, 7, 7),
##         byrow = TRUE, nrow = 4, ncol = 6
##         ), 
## )
## ##layout.show(7)
## plotit(dw.ud, ylim = c(4, 5.5)) # 1
## plotit(dw.uu, ylim = c(0.5, 3)) # 3
## plotit(dw.Du, ylim = c(0.5, 3)) # 2
## plotit(dw.Dd, ylim = c(3.5, 8.5)) # 4
## plotit(lat.Dd, ylim = c(4, 6.5)) # 5
## bracket <- c(-0.2, 0.2)
## plotit(lat.Dd, ylim = 6.1 + bracket, xlim = 0.91 + 2*0.1*bracket) # 6
## plotit(lat.Dd, ylim = 4.5 + bracket, xlim = 0.14 + 2*0.1*bracket) # 7
## if(save_plots) dev.off()

#### old code below here ####

## ht <- 15
## wd <- 15
## if(save_plots) {
##     pdf("./img/example-dw.pdf", height = ht, width = wd)
## } else { 
##     dev.new(height = ht, width = wd)
## }
## par(mfrow = c(2, 2), mai = c(1, 1, 0.1, 0.1))
## with(
##     list(
##         X = dw.Dd[[1]],
##         D = attr(dw.Dd, "bparam.vals")
##     ), {
##         mins <- apply(X, 1, min)
##         far1 <- 10
##         far2 <- 35
##         nearest <- which(mins < 3)[1] - 1
##         near <- nearest - 3
##                                         # marking lines; will need the SSDs for these lines
##         markD <- D[c(far1, far2, near, nearest)]
##         markX <- X[c(far1, far2, near, nearest), ]
##         print(apply(markX, 1, sd))
##                                         # plot        
##         matplot(
##             D, X, type = "l", lty = 1, lwd = 0.5, col = 9,
##             ylim = c(3.5, 8.5), xlab = "D", ylab = expression(x[i]), axes = FALSE, font.lab = 3, cex.lab = 2
##         )
##         axis(1, cex.axis = 2)
##         axis(2, cex.axis = 2)
##         segments(
##             x0 = markD,
##             y0 = apply(markX, 1, min) - 0.1,
##             y1 = apply(markX, 1, max) + 0.1,
##             col = 3, lty = 1, lwd = 4
##         )
##     }
## )
## with(
##     list(
##         X = dw.ud[[1]],
##         u = attr(dw.ud, "bparam.vals")
##     ), {
##         mins <- apply(X, 1, min)
##         far1 <- 10
##         far2 <- 35
##         nearest <- which(mins < 3)[1] - 1
##         near <- nearest - 3
##                                         # marking lines; will need the SSDs for these lines
##         marku <- u[c(far1, far2, near, nearest)]
##         markX <- X[c(far1, far2, near, nearest), ]
##         print(apply(markX, 1, sd))
##                                         # plot
##         matplot(
##             u, X, type = "l", lty = 1, lwd = 0.5, col = 9,
##             ylim = c(4, 5.5), xlab = "u", ylab = "", axes = FALSE, font.lab = 3, cex.lab = 2
##         )
##         axis(1, cex.axis = 2)
##         axis(2, cex.axis = 2)
##         segments(
##             x0 = marku,
##             y0 = apply(markX, 1, min) - 0.05,
##             y1 = apply(markX, 1, max) + 0.05,
##             col = 3, lty = 1, lwd = 4
##         )
##     }
## )
## with(
##     list(
##         X = dw.Du[[1]],
##         D = attr(dw.Du, "bparam.vals")
##     ), {
##         maxs <- apply(X, 1, max)
##         far1 <- 10
##         far2 <- 35
##         nearest <- which(maxs < 3)
##         nearest <- nearest[length(nearest)]
##         near <- nearest - 3
##                                         # marking lines; will need the SSDs for these lines
##         markD <- D[c(far1, far2, near, nearest)]
##         markX <- X[c(far1, far2, near, nearest), ]
##         print(apply(markX, 1, sd))
##                                         # plot        
##         matplot(
##             D, X, type = "l", lty = 1, lwd = 0.5, col = 9,
##             ylim = c(0.5, 3), xlim = c(0, 0.15),
##             xlab = "D", ylab = expression(x[i]), axes = FALSE, font.lab = 3, cex.lab = 2
##         )
##         axis(1, cex.axis = 2)
##         axis(2, cex.axis = 2)
##         segments(
##             x0 = markD,
##             y0 = apply(markX, 1, min) - 0.1,
##             y1 = apply(markX, 1, max) + 0.1,
##             col = 3, lty = 1, lwd = 4
##         )
##     }
## )
## with(
##     list(
##         X = dw.uu[[1]],
##         u = attr(dw.uu, "bparam.vals")
##     ), {
##         maxs <- apply(X, 1, max)
##         far1 <- 10
##         far2 <- 35
##         nearest <- which(maxs < 3)
##         nearest <- nearest[length(nearest)]
##         near <- nearest - 3
##                                         # marking lines; will need the SSDs for these lines
##         marku <- u[c(far1, far2, near, nearest)]
##         markX <- X[c(far1, far2, near, nearest), ]
##         print(apply(markX, 1, sd))
##                                         # plot
##         matplot(
##             u, X, type = "l", lty = 1, lwd = 0.5, col = 9,
##             ylim = c(0.5, 3), xlab = "u", ylab = "", axes = FALSE, font.lab = 3, cex.lab = 2
##         )
##         axis(1, cex.axis = 2)
##         axis(2, cex.axis = 2)
##         segments(
##             x0 = marku,
##             y0 = apply(markX, 1, min) - 0.05,
##             y1 = apply(markX, 1, max) + 0.05,
##             col = 3, lty = 1, lwd = 4
##         )
##     }
## )
## if(save_plots) dev.off()

## ht <- 7.5
## wd <- ht*3
## if(save_plots) {
##     pdf("./img/example-lattice-dw-D.pdf", height = ht, width = wd)
## } else {
##     dev.new(height = ht, width = wd)
## }
## test <- readRDS("./data/sims/doublewell-lattice-D-sde-50-down--5-1-TRUE.rds")
## par(mfrow = c(1, 3), mai = c(1, 1, 0.1, 0.1))
## matplot(attr(test, "bparam.vals"), test[[1]], type = "l", lty = 1, col = 9, lwd = 0.5,
##         ylim = c(4, 6.5), xlab = "D", ylab = expression(x[i]),
##         axes = FALSE, cex.lab = 2)
## axis(1, cex.axis = 2)
## axis(2, cex.axis = 2)
## matplot(attr(test, "bparam.vals"), test[[1]], type = "l", lty = 1, col = 9, lwd = 0.5,
##         xlim = c(0.12, 0.15), ylim = c(4.3, 4.7), xlab = "D", ylab = expression(x[i]),
##         axes = FALSE, cex.lab = 2)
## axis(1, cex.axis = 2)
## axis(2, cex.axis = 2)
## matplot(attr(test, "bparam.vals"),  test[[1]], type = "l", lty = 1, col = 9, lwd = 0.5,
##         xlim = c(0.88, 0.91), ylim = c(5.85, 6.25), xlab = "D", ylab = expression(x[i]),
##         axes = FALSE, cex.lab = 2)
## axis(1, cex.axis = 2)
## axis(2, cex.axis = 2)
## if(save_plots) dev.off()
