## New plan.
## Make as two separate figures here, then combine in Inkscape.
## Top row is like "example-method": the xi, then Moran's I, then SSD.
## Bottom row is the network figures and the colorbar.
## Use the xi plot from this file, with the four sample markers. 


library(sfsmisc)
library(igraph)
source("calc-functions.R")

save_plots <- TRUE # FALSE
palette("Okabe-Ito")

get_colors <- function(xi,
                       hcl.pal = NULL, n = NULL,
                       col.min = NULL, col.max = NULL,
                       rescale.cols = FALSE, xi.max = NULL, ...) {
    if(is.null(hcl.pal)) {
        cr <- colorRamp(c(col.min, col.max))
    } else {
        cr <- colorRamp(hcl.colors(n, hcl.pal, ...))
    }
    if(rescale.cols) {
        if(is.null(xi.max)) {
            colors <- cr(xi/max(xi))
        } else {
            colors <- cr(xi/xi.max)
        }
    } else {
        colors <- cr(xi) # all xi for SIS are between zero and one
    }
    cols <- apply(colors, 1, function(row) rgb(row[1], row[2], row[3], maxColorValue = 255))
    cols
}

## draw_spatial_correlation <- function(g, xi, col.min = "white", col.max = 6, ...) {
draw_spatial_correlation <- function(g, xi, ...) {
    plot(g, vertex.size = 10, vertex.label = "",
         vertex.frame.color = adjustcolor(9, 0.5), 
         ##vertex.color = get_colors(xi, col.min, col.max, rescale.cols = FALSE),
         vertex.color = get_colors(xi, hcl.pal = "Spectral", n = 5, rescale.cols = TRUE, xi.max = NULL, rev = TRUE),
         edge.color = adjustcolor(9, 0.5), ...)
}

networks <- readRDS("./data/networks.rds")

SIS.Du <- readRDS("./data/sims/SIS-drug-D-sde-50-up-0-1-TRUE.rds")
SIS.Dd <- readRDS("./data/sims/SIS-drug-D-sde-50-down-0-1-TRUE.rds")
dw.Dd <- readRDS("./data/sims/doublewell-drug-D-sde-50-down--5-1-TRUE.rds")
## try showing catlins and canton
canton <- readRDS("./data/sims/doublewell-canton-D-sde-50-down--5-1-TRUE.rds")
catlins <- readRDS("./data/sims/doublewell-catlins-D-sde-50-down--5-1-TRUE.rds")

sim <- SIS.Du # catlins canton
g <- networks[["drug"]] # catlins canton
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)
X <- sim[[1]]
xi.max <- max(X)
bparam.vals <- attr(sim, "bparam.vals")
direction <- attr(sim, "direction")
model <- attr(sim, "model")
bpoint <- switch(model, doublewell = 3, SIS = 0.01)

edges <- switch(direction, down = apply(X, 1, min), up = apply(X, 1, max))
far1 <- 10
far2 <- 35
nearest <- switch(direction, down = which(edges < bpoint)[1] - 1, up = which(edges > bpoint)[1] - 1)
near <- nearest - 3

                                        # data for marking lines, as matrix
markvec <- c(far1, far2, near, nearest)
mark.bpv <- bparam.vals[markvec]
mark.X <- X[markvec, ]
pad <- 0.025*rowMeans(mark.X)
xi.max <- max(mark.X)

pht <- pwd <- 4

wd <- 3.5*pwd
ht <- 1.5*pwd
if(save_plots) {
    pdf("./img/example-SIS-up-top.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
layoutmatrix <- matrix(1:3, byrow = TRUE, nrow = 1)
thelayout <- layout(layoutmatrix, widths = rep(lcm(2.25*pwd), 3), heights = lcm(2.25*pht))

df <- promote_df(sim)
A <- attr(df, "params")$A
D <- attr(df, "bparam.vals")

idx <- get_idx(df)
xlim <- c(0, 0.15)#range(bparam.vals[idx])

                                        # The xi, panel (a)
par(mar = c(5, 5, 0.5, 0.5))
matplot(bparam.vals, X, type = "l", lty = 1, col = 9,
        xlab = "D", ylab = expression(x[i]), axes = FALSE, cex.lab = 2, xlim = xlim, ylim = c(0, 1e-2))
eaxis(1, cex.axis = 2)
eaxis(2, cex.axis = 2)
segments(
    x0 = mark.bpv, y0 = apply(mark.X, 1, min) - pad, y1 = apply(mark.X, 1, max) + pad,
        col = 2, lty = 1, lwd = 4
)


moranI <- apply(df, 1, global_moran, A = A)[idx]
ssd <- apply(df, 1, sd)[idx]

print(tau.I <- get_tau(df, moranI, TRUE))
print(tau.ssd <- get_tau(df, ssd, TRUE))

n <- 5
farsample <- get_samples(df, "far", n = n)
nearsample <- get_samples(df, "near", n = n)
farmodel.I <- get_slope(df, moranI, which = "far", n = n, return.model = TRUE)
nearmodel.I <- get_slope(df, moranI, which = "near", n = n, return.model = TRUE)
farmodel.sd <- get_slope(df, ssd, which = "far", n = n, return.model = TRUE)
nearmodel.sd <- get_slope(df, ssd, which = "near", n = n, return.model = TRUE)

                                        # Moran's I, panel (b)
plot(D[idx], moranI[idx], type = "l", xlim = xlim, #xlim = range(D),
     lty = 1, lwd = 2, col = adjustcolor(2, 0.5),
     xlab = "", ylab = "", axes = FALSE)
lines(D[farsample], predict(farmodel.I), lwd = 4, col = 2, lty = 1)
lines(D[nearsample], predict(nearmodel.I), lwd = 4, col = 2, lty = 1)
eaxis(1, cex.axis = 1.75)
eaxis(2, cex.axis = 1.75)
title(xlab = attr(df, "bparam"), cex.lab = 2)
##title(ylab = "Moran's I", cex.lab = 2, line = 3.2)
legend("bottomright", col = 2, lwd = 2, lty = 1, legend = "Moran's I", bty = "n", cex = 1.75)

                                        # SSD, panel (c)
plot(D[idx], ssd[idx], type = "l", xlim = xlim, #xlim = range(D),
     lty = 1, lwd = 2, col = adjustcolor(3, 0.5),
     xlab = "", ylab = "", axes = FALSE)
lines(D[farsample], predict(farmodel.sd), lwd = 4, col = 3, lty = 1)
lines(D[nearsample], predict(nearmodel.sd), lwd = 4, col = 3, lty = 1)
eaxis(1, cex.axis = 1.75)
eaxis(2, cex.axis = 1.75)
title(xlab = attr(df, "bparam"), cex.lab = 2)
##title(ylab = "Spatial standard deviation", cex.lab = 2, line = 3.2)
legend("bottomright", col = 3, lwd = 2, lty = 1, legend = "Std. dev.", bty = "n", cex = 1.75)

if(save_plots) dev.off()

### new figure


wd <- 4.5*pwd#3.25*pwd
ht <- 1.5*pwd#2*pht
if(save_plots) {
    pdf("./img/example-SIS-up-bottom.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
layoutmatrix <- matrix(1:4, nrow = 1, byrow = TRUE)
thelayout <- layout(layoutmatrix, widths = rep(lcm(2.25*pwd), 4), heights = lcm(2.25*pht))
## layout.show(thelayout)

line1 <- -1.1
line2 <- -2.2
line3 <- -3.3
set.seed(123)
lyt <- layout_nicely(g)
par(mar = rep(0, 4))

draw_spatial_correlation(g, X[far1, ], layout = lyt)
mtext("Farthest", line = line1, adj = 0.02, col = 1)
mtext(paste(round(global_moran(X[far1, ], A), 3)), line = line2, adj = 0.02, col = 2)
mtext(paste(round(sd(X[far1, ]), 3)), line = line3, adj = 0.02, col = 3)

draw_spatial_correlation(g, X[far2, ], layout = lyt)
mtext("Far", line = line1, adj = 0.02, col = 1)
mtext(paste(round(global_moran(X[far2, ], A), 3)), line = line2, adj = 0.02, col = 2)
mtext(paste(round(sd(X[far2, ]), 3)), line = line3, adj = 0.02, col = 3)

draw_spatial_correlation(g, X[near, ], layout = lyt)
mtext("Near", line = line1, adj = 0.02, col = 1)
mtext(paste(round(global_moran(X[near, ], A), 3)), line = line2, adj = 0.02, col = 2)
mtext(paste(round(sd(X[near, ]), 3)), line = line3, adj = 0.02, col = 3)

draw_spatial_correlation(g, X[nearest, ], layout = lyt)
mtext("Nearest", line = line1, adj = 0.02, col = 1)
mtext(paste(round(global_moran(X[nearest, ], A), 3)), line = line2, adj = 0.02, col = 2)
mtext(paste(round(sd(X[nearest, ]), 3)), line = line3, adj = 0.02, col = 3)

if(save_plots) dev.off()

### color bar
colorbar <- function(colors, rng, title = "", ...) {
    ##crp <- colorRampPalette(colors)
    z <- matrix(1:100, nrow = 1)
    x <- 1
    y <- seq(rng[1], rng[2], length.out = 100)
                                        ##par(mar = c(3, 0, 6, 3.5) + .5)
    par(mar = c(0, 12, 0, 12))
    image(x, y, z, col = colors,#crp(100),
          axes = FALSE, xlab = "", ylab = "")
    ##eaxis(4, at = pretty(rng), labels = pretty(rng), ...)
    eaxis(2, ...)
    mtext(title, 4, line = 2.5)
    box()
}

wd <- 4.5*pwd#3.25*pwd
ht <- 1.5*pwd#2*pht
if(save_plots) {
    pdf("./img/example-SIS-up-bottom-colorbars.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
layoutmatrix <- matrix(1:4, nrow = 1, byrow = TRUE)
thelayout <- layout(layoutmatrix, widths = rep(lcm(2.25*pwd), 4), heights = lcm(2.25*pht))
## layout.show(thelayout)

colorbar(hcl.colors(100, "Spectral", rev = TRUE), range(X[far1, ]), cex.axis = 2)
colorbar(hcl.colors(100, "Spectral", rev = TRUE), range(X[far2, ]), cex.axis = 2)
colorbar(hcl.colors(100, "Spectral", rev = TRUE), range(X[near, ]), cex.axis = 2)
colorbar(hcl.colors(100, "Spectral", rev = TRUE), range(X[nearest, ]), cex.axis = 2)

if(save_plots) dev.off()

## ht <- 15/2
## wd <- 15
## if(save_plots) {
##     pdf("./img/example-SIS.pdf", height = ht, width = wd)
## } else { 
##     dev.new(height = ht, width = wd)
## }
## par(mfrow = c(1, 2), mai = c(1, 1, 0.1, 0.1))
## with(
##     list(
##         X = SIS.Du[[1]],
##         D = attr(SIS.Du, "bparam.vals")
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
##         print(apply(markX, 1, global_moran, A))
##                                         # plot
##         matplot(
##             D, X, type = "l", lty = 1, lwd = 0.5, col = 9,
##             ylim = c(0, 0.01), xlim = c(0, 0.2),
##             xlab = "D", ylab = expression(x[i]), axes = FALSE, font.lab = 3, cex.lab = 2
##         )
##         axis(1, cex.axis = 2)
##         axis(2, cex.axis = 2)
##         segments(
##             x0 = markD,
##             y0 = apply(markX, 1, min) - 0.005,
##             y1 = apply(markX, 1, max) + 0.005,
##             col = 3, lty = 1, lwd = 4
##         )
##     }
## )
