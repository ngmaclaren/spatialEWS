library(sfsmisc)
library(igraph)
source("calc-functions.R")
save_plots <- FALSE # TRUE
palette("Okabe-Ito")

dw.Dd <- readRDS("./data/sims/doublewell-drug-D-sde-50-down--5-1-TRUE.rds")
dw.Du <- readRDS("./data/sims/doublewell-drug-D-sde-50-up-0-1-TRUE.rds")

g <- readRDS("./data/networks.rds")[["drug"]]
## write_graph(g, "./img/drug.graphml", "graphml")

colors <- list(
    base = 1,
    moran = 2,
    cv = 3, #ssd = 3,
    skew = 4,
    kurt = 5,
    xi = 9
)

plotit <- function(sim, ...) {
    df <- promote_df(sim)
    A <- attr(df, "params")$A
    D <- attr(df, "bparam.vals")

    idx <- get_idx(df)
    
    moranI <- apply(df, 1, global_moran, A = A)[idx]
    cv <- apply(df, 1, CV)[idx] # ssd <- apply(df, 1, sd)[idx]
    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew
    kurt <- apply(df, 1, moments::kurtosis)

    taus <- list(
        I_M = get_tau(df, moranI, TRUE),
        cv = get_tau(df, cv, TRUE), # ssd = get_tau(df, ssd, TRUE),
        skew = get_tau(df, skew, TRUE),
        kurt = get_tau(df, kurt, TRUE)
    )
    print(taus)

    print(classify(df, moranI, n = 5))
    print(classify(df, cv, n = 5)) # print(classify(df, ssd, n = 5))
    print(classify(df, skew, n = 5))
    print(classify(df, kurt, n = 5))

    pos <- 5.25
    lwd <- 3

    par(mar = c(5, 5, 0.5, 20)) # lots of space for additional axes
    matplot(D, df, type = "l", lty = 1, lwd = 0.5, col = colors$xi, xlab = "", ylab = "", axes = FALSE)
    eaxis(1, cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    title(xlab = attr(df, "bparam"), cex.lab = 2)
    title(ylab = expression(x[i]), cex.lab = 2)#, line = 3)

    par(new = TRUE)
    plot(D[idx], moranI[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$moran,
         xlab = "", ylab = "", axes = FALSE)
    eaxis(
        4, cex.axis = 1.75, col = colors$moran, col.axis = colors$moran,
        small.args = list(col = colors$moran)
    )

    par(new = TRUE)
    plot(D[idx], cv[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$cv, #ssd,
         xlab = "", ylab = "", axes = FALSE)
    eaxis(
        4, cex.axis = 1.75, col = colors$cv, line = pos, col.axis = colors$cv, # ssd
        small.args = list(col = colors$cv, line = pos)
    )

    par(new = TRUE)
    plot(D[idx], skew[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$skew,
         xlab = "", ylab = "", axes = FALSE)
    eaxis(
        4, cex.axis = 1.75, col = colors$skew, line = pos*2, col.axis = colors$skew,
        small.args = list(col = colors$skew, line = pos*2)
    )

    par(new = TRUE)
    plot(D[idx], kurt[idx], type = "l", xlim = range(D), lty = 1, lwd = lwd, col = colors$kurt,
         xlab = "", ylab = "", axes = FALSE)
    eaxis(
        4, cex.axis = 1.75, col = colors$kurt, line = pos*3, col.axis = colors$kurt,
        small.args = list(col = colors$kurt, line = pos*3)
    )
}

ht <- 10
wd <- 9
if(save_plots) {
    pdf("./img/example-method.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(2, 1))

## plot(g, vertex.label = "", vertex.size = 5, vertex.color = 5, vertex.frame.color = NA,
##      edge.color = 8)

plotit(dw.Du)
plotit(dw.Dd)

if(save_plots) dev.off()
