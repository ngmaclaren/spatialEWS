library(sfsmisc)
library(igraph)
source("calc-functions.R")
save_plots <- TRUE # FALSE
palette("Okabe-Ito")

colors <- list(
    base = 1,
    moran = 2,
    cv = 3, #ssd = 3,
    skew = 4,
    kurt = 8,
    xi = 9
)

plotit <- function(sim, ...) {
    df <- promote_df(sim)
    A <- attr(df, "params")$A
    D <- attr(df, "bparam.vals")

    idx <- get_idx(df)
    
    moranI <- apply(df, 1, global_moran, A = A)[idx]
    cv <- apply(df, 1, CV)[idx]
    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew
    kurt <- apply(df, 1, moments::kurtosis)

    taus <- list(
        I_M = get_tau(df, moranI, TRUE),
        cv = get_tau(df, cv, TRUE),
        skew = get_tau(df, skew, TRUE),
        kurt = get_tau(df, kurt, TRUE)
    )
    classifications <- list(
        I_M = classify(df, moranI, n = 5),
        cv = classify(df, cv, n = 5),
        skew = classify(df, skew, n = 5),
        kurt = classify(df, kurt, n = 5)
    )

    labelsize <- 2.25
    lwd <- 2
    
    layoutmatrix <- matrix(c(
        1, 1, 2, 3,
        1, 1, 4, 5), byrow = TRUE, nrow = 2)
    layout(layoutmatrix)

    par(mar = c(5, 5, 1, 1), pty = "s")
    matplot(D, df, type = "l", lty = 1, lwd = 0.5, col = colors$xi, xlab = "", ylab = "", axes = FALSE)
    box()
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(xlab = attr(df, "bparam"), cex.lab = 1.15*labelsize)
    title(ylab = expression(x[i]), cex.lab = 1.15*labelsize)

    plot(D[idx], moranI[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$moran,
         xlab = "", ylab = "", axes = FALSE)
    mtext(
        paste0("tau = ", round(taus$I_M, 2), "\n", tosentence(classifications$I_M)),
        cex = 0.75, col = colors$moran, line = -3, adj = 0.98
    )
    box()
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(xlab = attr(df, "bparam"), cex.lab = 1.15*labelsize)
    ##title(ylab = "Moran's I", cex.lab = 1.15*labelsize)

    plot(D[idx], cv[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$cv, #ssd,
         xlab = "", ylab = "", axes = FALSE)
    mtext(
        paste0("tau = ", round(taus$cv, 2), "\n", tosentence(classifications$cv)),
        cex = 0.75, col = colors$cv, line = -3, adj = 0.98
    )
    box()
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(xlab = attr(df, "bparam"), cex.lab = 1.15*labelsize)
    ##title(ylab = "CV", cex.lab = 1.15*labelsize)

    plot(D[idx], skew[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$skew,
         xlab = "", ylab = "", axes = FALSE)
    mtext(
        paste0("tau = ", round(taus$skew, 2), "\n", tosentence(classifications$skew)),
        cex = 0.75, col = colors$skew, line = -3, adj = 0.98
    )
    box()
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(xlab = attr(df, "bparam"), cex.lab = 1.15*labelsize)
    ##title(ylab = "Skewness", cex.lab = 1.15*labelsize)

    plot(D[idx], kurt[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$kurt,
         xlab = "", ylab = "", axes = FALSE)
    mtext(
        paste0("tau = ", round(taus$kurt, 2), "\n", tosentence(classifications$kurt)),
        cex = 0.75, col = colors$kurt, line = -3, adj = 0.98
    )
    box()
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(xlab = attr(df, "bparam"), cex.lab = 1.15*labelsize)

}

dw.Dd <- readRDS("./data/sims/doublewell-drug-D-sde-50-down--5-1-TRUE.rds")
dw.Du <- readRDS("./data/sims/doublewell-drug-D-sde-50-up-0-1-TRUE.rds")

g <- readRDS("./data/networks.rds")[["drug"]]

ht <- 5
wd <- 10
if(save_plots) {
    pdf("./img/example-method-top.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
plotit(dw.Du)
if(save_plots) dev.off()

if(save_plots) {
    pdf("./img/example-method-bottom.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
plotit(dw.Dd)
if(save_plots) dev.off()

sis.Dd <- readRDS("./data/sims/SIS-lattice-D-sde-50-down-0-1-TRUE.rds")
sis.Du <- readRDS("./data/sims/SIS-lattice-D-sde-50-up-0-1-TRUE.rds")

g <- readRDS("./data/networks.rds")[["lattice"]]

if(save_plots) {
    pdf("./img/example-lattice-top.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
plotit(sis.Du)
if(save_plots) dev.off()

if(save_plots) {
    pdf("./img/example-lattice-bottom.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
plotit(sis.Dd)
if(save_plots) dev.off()
