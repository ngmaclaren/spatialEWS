library(sfsmisc)
library(igraph)
source("calc-functions.R")
save_plots <- FALSE # TRUE
palette("Okabe-Ito")

dw.Dd <- readRDS("./data/sims/doublewell-drug-D-sde-50-down--5-1-TRUE.rds")
dw.Du <- readRDS("./data/sims/doublewell-drug-D-sde-50-up-0-1-TRUE.rds")

g <- readRDS("./data/networks.rds")[["drug"]]
## write_graph(g, "./img/drug.graphml", "graphml")

plotit <- function(sim, n = 3, ...) {
    df <- promote_df(sim)
    A <- attr(df, "params")$A
    D <- attr(df, "bparam.vals")

    idx <- get_idx(df)
    
    moranI <- apply(df, 1, global_moran, A = A)[idx]
    ssd <- apply(df, 1, sd)[idx]

    print(tau.I <- get_tau(df, moranI, TRUE))
    print(tau.ssd <- get_tau(df, ssd, TRUE))

    farsample <- get_samples(df, "far", n = n)
    nearsample <- get_samples(df, "near", n = n)
    farmodel.I <- get_slope(df, moranI, which = "far", n = n, return.model = TRUE)
    nearmodel.I <- get_slope(df, moranI, which = "near", n = n, return.model = TRUE)
    farmodel.sd <- get_slope(df, ssd, which = "far", n = n, return.model = TRUE)
    nearmodel.sd <- get_slope(df, ssd, which = "near", n = n, return.model = TRUE)

    par(mar = c(6, 6, 0.5, 0.5))
    matplot(D, df, type = "l", lty = 1, lwd = 0.5, col = 9, xlab = "", ylab = "", axes = FALSE)
    eaxis(1, cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    title(xlab = attr(df, "bparam"), cex.lab = 2)
    title(ylab = expression(x[i]), cex.lab = 2, line = 3.2)

    plot(D[idx], moranI[idx], type = "l", xlim = range(D), lty = 1, lwd = 2, col = adjustcolor(2, 0.5),
         xlab = "", ylab = "", axes = FALSE)
    lines(D[farsample], predict(farmodel.I), lwd = 4, col = 2, lty = 1)
    lines(D[nearsample], predict(nearmodel.I), lwd = 4, col = 2, lty = 1)
    eaxis(1, cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    title(xlab = attr(df, "bparam"), cex.lab = 2)
    ##title(ylab = "Moran's I", cex.lab = 2, line = 3.2)
    legend("bottomright", col = 2, lwd = 2, lty = 1, legend = "Moran's I", bty = "n", cex = 1.75)

    plot(D[idx], ssd[idx], type = "l", xlim = range(D), lty = 1, lwd = 2, col = adjustcolor(3, 0.5),
         xlab = "", ylab = "", axes = FALSE)
    lines(D[farsample], predict(farmodel.sd), lwd = 4, col = 3, lty = 1)
    lines(D[nearsample], predict(nearmodel.sd), lwd = 4, col = 3, lty = 1)
    eaxis(1, cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    title(xlab = attr(df, "bparam"), cex.lab = 2)
    ##title(ylab = "Spatial standard deviation", cex.lab = 2, line = 3.2)
    legend("bottomright", col = 3, lwd = 2, lty = 1, legend = "Std. dev.", bty = "n", cex = 1.75)
}

ht <- 8
wd <- 12
if(save_plots) {
    pdf("./img/example-method.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
layoutmatrix <- matrix(
    c(1, 2, 3,
      4, 5, 6),
    byrow = TRUE, nrow = 2
)
layout(layoutmatrix)
## layout.show(6)

## plot(g, vertex.label = "", vertex.size = 5, vertex.color = 5, vertex.frame.color = NA,
##      edge.color = 8)

plotit(dw.Dd, n = 5)
plotit(dw.Du, n = 5)

if(save_plots) dev.off()
