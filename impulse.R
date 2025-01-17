library(igraph)
library(sdn)
library(moments)
library(sfsmisc)
source("calc-functions.R")
palette("Okabe-Ito")

set.seed(1)
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
A <- as_adj(g, "both", sparse = FALSE)
AL <- as_adj_list(g, "all")
N <- vcount(g)

model <- doublewell
params <- c(.doublewell, list(AL = AL))
deltaT <- 0.01

simtime.pre <- 50 
simtime.post <- 51

xinit <- rep(params$xinit.high, N)
control <- list(times = 0:simtime.pre, deltaT = deltaT)

params$u <- -2
system.time(X.pre <- sde(xinit, control$times, model, params, control))

control$times <- 0:simtime.post
params$u <- -4
system.time(X.post <- sde(X.pre[nrow(X.pre), -1], control$times, model, params, control))

X.post[, "time"] <- X.post[, "time"] + simtime.pre

X <- rbind(X.pre, X.post)
attr(X, "direction") <- "down"
attr(X, "dynamics") <- "doublewell"
attr(X, "deltaT") <- deltaT
basin <- basins[["doublewell"]]

##idx <- which(X[, "time"] >= 5 & X[, "time"] <= simtime.pre+simtime.post-1)
##mark <- which(X[, "time"] == simtime.pre)
idx <- get_idx_time(X)
mark <- max(idx)
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
    pdf(paste0("./img/impulse-", network, ".pdf"), height = ht, width = wd)
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
    X[plotidx, 1], X[plotidx, -1], type = "l", col = 1, lty = 1, lwd = 0.5,
    xlab = "", ylab = "", axes = FALSE, xaxs = "i",
    ylim = c(0, ceiling(max(X[, -1]))), yaxs = "i"
)
## abline(v = X[mark, "time"], lty = 2, lwd = 2, col = colors$mark) # maybe?
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
title(ylab = "x", cex.lab = 1.15*labelsize, font.lab = 3)
                                        # 2: Moran's I
plot(X[idx, "time"], moranI[idx], type = "l", col = colors$moran, lty = 1, lwd = 2,
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", xlim = range(X[, "time"]))
##abline(v = X[mark, "time"], lty = 3, lwd = 2, col = colors$mark) # simtime.pre
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
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", xlim = range(X[, "time"]))
## abline(v = X[mark, "time"], lty = 3, lwd = 2, col = colors$mark)
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
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", xlim = range(X[, "time"]))
## abline(v = X[mark, "time"], lty = 3, lwd = 2, col = colors$mark)
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
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", xlim = range(X[, "time"]))
##abline(v = X[mark, "time"], lty = 3, lwd = 2, col = colors$mark)
mtext(
    paste0("tau = ", round(taus$kurt, 2), "\n", tosentence(classifications$kurt)),
    cex = 0.75, col = colors$kurt, line = -3, adj = 0.98
)
box()
eaxis(1, cex.axis = labelsize)
eaxis(2, cex.axis = labelsize)
title(xlab = "t", cex.lab = 1.15*labelsize, font.lab = 3)
## title(ylab = "Kurtosis", cex.lab = 1.15*labelsize)
if(save_plots) dev.off()
