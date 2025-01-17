library(igraph)
library(sdn)
library(moments)
library(sfsmisc)
source("calc-functions.R")
palette("Okabe-Ito")

save_plots <- TRUE # FALSE
networks <- readRDS("./data/networks.rds")
colors <- list(
    base = 1,
    moran = 2,
    cv = 3,
    skew = 4,
    kurt = 8,
    xi = 1,#9, # to differentiate from equilibrium values?
    mark = 7
)

## network <- "catlins"
network <- "ug_village"
bparam <- "u"
direction <- "down"

g <- networks[[network]]
A <- as_adjacency_matrix(g, "both", sparse = FALSE)
AL <- as_adj_list(g, "all")
N <- vcount(g)

model <- doublewell
params <- c(.doublewell, list(AL = AL))
deltaT <- 0.01

simtime <- 50

xinit <- rep(params$xinit.high, N)
control <- list(times = 0:simtime, deltaT = deltaT)

if(network == "catlins") set.seed(3) else if(network == "ug_village") set.seed(4)
params$u <- -3
if(network == "catlins") params$sigma <- 0.25 else if(network == "ug_village") params$sigma <- 0.2
system.time(X <- sde(xinit, control$times, model, params, control))
attr(X, "direction") <- "down"
attr(X, "dynamics") <- "doublewell"
attr(X, "deltaT") <- deltaT
basin <- basins[["doublewell"]]

idx <- get_idx_time(X)
mark <- max(idx)
zoom <- X[c(mark - 500, mark + 500), "time"]
unzoom <- range(X[, "time"])
moranI <- apply(X[, -1], 1, global_moran, A = A)
cv <- apply(X[, -1], 1, CV)
skew <- -apply(X[, -1], 1, skewness) # negate b/c descending sims
kurt <- apply(X[, -1], 1, kurtosis)

taus <- list(
    I_M = get_tau_time(X, moranI),
    cv = get_tau_time(X, cv),
    skew = get_tau_time(X, skew),
    kurt = get_tau_time(X, kurt)
)

classifications <- list(
    I_M = classify_time(X, moranI),
    cv = classify_time(X, cv),
    skew = classify_time(X, skew),
    kurt = classify_time(X, kurt)
)

ht <- 5
wd <- 10
labelsize <- 2.25
if(save_plots) {
    pdf(paste0("./img/stochastic-", network, ".pdf"), height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
layoutmatrix <- matrix(c(
    1, 1, 2, 3,
    1, 1, 4, 5), byrow = TRUE, nrow = 2)
layout(layoutmatrix)
## layout.show(5)
par(mar = c(5, 5, 1, 1), pty = "s")
                                        # 1
plotidx <- seq(1, nrow(X), by = 1/deltaT)#(1/deltaT)/10)
matplot(
    X[plotidx, 1], X[plotidx, -1], type = "l", col = colors$xi, lty = 1, lwd = 0.5,
    xlab = "", ylab = "", axes = FALSE, xaxs = "i",
    ylim = c(0, ceiling(max(X[, -1]))), yaxs = "i"
)
abline(v = X[mark, "time"], lty = 2, lwd = 2, col = colors$mark)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
title(ylab = "x", cex.lab = 1.15*labelsize, font.lab = 3)
                                        # 2: Moran's I
plot(X[idx, "time"], moranI[idx], type = "l", col = colors$moran, lty = 1, lwd = 2,
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", xlim = unzoom)
abline(v = X[mark, "time"], lty = 2, lwd = 2, col = colors$mark)
mtext(
    paste0("tau = ", round(taus$I_M, 2), "\n", tosentence(classifications$I_M)),
    cex = 0.75, col = colors$moran, line = -3, adj = 0.98
)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
##title(ylab = "Moran's I", cex.lab = 1.15*labelsize)
                                        # 3: CV
plot(X[idx, "time"], cv[idx], type = "l", col = colors$cv, lty = 1, lwd = 2,
     xlab = "", ylab = "", axes = FALSE, xaxs = "i",
     xlim = unzoom)
abline(v = X[mark, "time"], lty = 2, lwd = 2, col = colors$mark)
mtext(
    paste0("tau = ", round(taus$cv, 2), "\n", tosentence(classifications$cv)),
    cex = 0.75, col = colors$cv, line = -3, adj = 0.98
)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
##title(ylab = "CV", cex.lab = 1.15*labelsize)
                                        # 4: g_1'
plot(X[idx, "time"], skew[idx], type = "l", col = colors$skew, lty = 1, lwd = 2,
     xlab = "", ylab = "", axes = FALSE, xaxs = "i",
     xlim = unzoom)
abline(v = X[mark, "time"], lty = 2, lwd = 2, col = colors$mark)
mtext(
    paste0("tau = ", round(taus$skew, 2), "\n", tosentence(classifications$skew)),
    cex = 0.75, col = colors$skew, line = -3, adj = 0.98
)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
##title(ylab = "Skewness", cex.lab = 1.15*labelsize)
                                        # 5: g_2
plot(X[idx, "time"], kurt[idx], type = "l", col = colors$kurt, lty = 1, lwd = 2,
     xlab = "", ylab = "", axes = FALSE, xaxs = "i",
     xlim = unzoom)
abline(v = X[mark, "time"], lty = 2, lwd = 2, col = colors$mark)
mtext(
    paste0("tau = ", round(taus$kurt, 2), "\n", tosentence(classifications$kurt)),
    cex = 0.75, col = colors$kurt, line = -3, adj = 0.98
)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
##title(ylab = "Kurtosis", cex.lab = 1.15*labelsize)
if(save_plots) dev.off()

#which(apply(X[, -1], 1, min) < 3)[1]
## cimarks <- which(X[, "time"] >= 5 & X[, "time"] <= 15)
## let's make a 95% confidence band around each EWS mean
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = 3,
##          lwd = 2, col = 9, lty = 1)
## abline(h = 3, col = 9, lwd = 0.5, lty = 2)
## bands <- list(
##     moranI = quantile(moranI[cimarks  - (cimarks[1] - 1)], probs = c(0.025, 0.975)),
##     cv = quantile(cv[cimarks  - (cimarks[1] - 1)], probs = c(0.025, 0.975)),
##     skew = quantile(skew[cimarks  - (cimarks[1] - 1)], probs = c(0.025, 0.975)),
##     kurt = quantile(kurt[cimarks  - (cimarks[1] - 1)], probs = c(0.025, 0.975))
## )
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$moranI[1],
##          lwd = 2, col = 9, lty = 1)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$moranI[2],
##          lwd = 2, col = 9, lty = 1)
## abline(h = bands$moranI, col = 9, lwd = 0.5, lty = 2)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$cv[1],
##          lwd = 2, col = 9, lty = 1)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$cv[2],
##          lwd = 2, col = 9, lty = 1)
## abline(h = bands$cv, col = 9, lwd = 0.5, lty = 2)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$skew[1],
##          lwd = 2, col = 9, lty = 1)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$skew[2],
##          lwd = 2, col = 9, lty = 1)
## abline(h = bands$skew, col = 9, lwd = 0.5, lty = 2)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$kurt[1],
##          lwd = 2, col = 9, lty = 1)
## segments(x0 = X[cimarks[1], "time"], x1 = X[cimarks[length(cimarks)], "time"], y0 = bands$kurt[2],
##          lwd = 2, col = 9, lty = 1)
## abline(h = bands$kurt, col = 9, lwd = 0.5, lty = 2)
