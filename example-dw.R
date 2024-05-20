## Use layout() to put all plots in the same file
## top row is three "good" plots
## bottom row is one bad plot, then lattice with two sections zoomed in
##
## For this and any other plots, replace the basic axes with sfsmisc axes


## choose a network: drug
## use stored sims

## Example #1: coupled double-well, showing where SSD works and where it doesn't.

## Example 2 should highlight when Moran's I does well, but I don't necessarily understand what's going on there
## One idea for Moran's I is to color the nodes by x_i. Do we see them visually becoming more patchy or similar as we approach the bifurcation?
library(sfsmisc)
save_plots <- TRUE # FALSE
palette("Okabe-Ito")

plotit <- function(sim, col.xi = 9, col.ssd = 3, pad.ssd = NULL, labelsize = NULL, ...) {
                                        # assumes palette is Okabe-Ito
    X <- sim[[1]]
    bparam.vals <- attr(sim, "bparam.vals")
    direction <- attr(sim, "direction")

                                        # Pick the bparam vals at which to compute ssd
    edges <- switch(
        direction,
        down = apply(X, 1, min),
        up = apply(X, 1, max)
    )
    far1 <- 10
    far2 <- 35
    nearest <- switch(
        direction,
        down = which(edges < 3)[1] - 1, # !! hard coding separatrix because in this file only considering dw
        up = which(edges > 3)[1] - 1
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

    print(apply(mark.X, 1, sd))
}

dw.Dd <- readRDS("./data/sims/doublewell-drug-D-sde-50-down--5-1-TRUE.rds")
dw.ud <- readRDS("./data/sims/doublewell-drug-u-sde-50-down-0-1-TRUE.rds")
dw.Du <- readRDS("./data/sims/doublewell-drug-D-sde-50-up-0-1-TRUE.rds")
dw.uu <- readRDS("./data/sims/doublewell-drug-u-sde-50-up-0-1-TRUE.rds")
lat.Dd <- readRDS("./data/sims/doublewell-lattice-D-sde-50-down--5-1-TRUE.rds")

pht <- pwd <- 5
wd <- 3*pwd
ht <- 2*pwd
if(save_plots) {
    pdf("./img/example-dw-v3.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
lyt <- layout(
    matrix(c(
        1, 1, 2, 2, 3, 3,
        1, 1, 2, 2, 3, 3,
        4, 4, 5, 5, 6, 6,
        4, 4, 5, 5, 7, 7),
        byrow = TRUE, nrow = 4, ncol = 6
        ), 
)
##layout.show(7)
plotit(dw.ud, ylim = c(4, 5.5)) # 1
plotit(dw.uu, ylim = c(0.5, 3)) # 3
plotit(dw.Du, ylim = c(0.5, 3)) # 2
plotit(dw.Dd, ylim = c(3.5, 8.5)) # 4
plotit(lat.Dd, ylim = c(4, 6.5)) # 5
bracket <- c(-0.2, 0.2)
plotit(lat.Dd, ylim = 6.1 + bracket, xlim = 0.91 + 2*0.1*bracket) # 6
plotit(lat.Dd, ylim = 4.5 + bracket, xlim = 0.14 + 2*0.1*bracket) # 7
if(save_plots) dev.off()

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
