## dolphin network
## D, down
## D = 0.4

## use sdn
library(igraph)
library(sdn)
library(runner)

save_plots <- TRUE # FALSE
palette("Okabe-Ito")

networks <- readRDS("./data/networks.rds")

g <- networks$dolphin
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)

xinit <- rep(.doublewell$xinit.high, N)
times <- 0:50
deltaT <- 0.01
control <- list(times = times, deltaT = deltaT)
params <- c(.doublewell, list(A = A))
params$u <- -5
## X <- sde(xinit, control$times, doublewell, params, control)
##summary(colMeans(tail(X, 20)))
## time_ev(X)

Ds <- c(0.6, 0.41, 0.377)
res <- lapply(Ds, function(D) {
    params$D <- D
    sde(xinit, control$times, doublewell, params, control)
})

examtimes <- seq(25, 50, by = deltaT)/deltaT

sds <- lapply(res, function(X) runner(X[examtimes, -1], sd, k = 1, na_pad = TRUE)) # this doesn't need runner()...
lim.sds <- range(unlist(sds))
lim.X <- range(unlist(lapply(res, function(X) as.numeric(X[, -1])))) 

if(save_plots) {
    pdf("./img/converge-diverge.pdf", height = 5, width = 17)
} else {
    dev.new(height = 5, width = 17)
}
par(mar = c(4, 4, 1, 4), mfcol = c(1, 3))
for(i in seq_along(res)) {
    time_ev(res[[i]], main = paste("D =", Ds[i]), ylim = lim.X)
    par(new = TRUE)
    plot(examtimes, sds[[i]], xlim = c(0, 50)/deltaT, ylim = lim.sds,
         type = "l", col = 3, axes = FALSE, xlab = "", ylab = "")
    axis(4, at = pretty(unlist(sds)), labels = pretty(unlist(sds)), col = 3, col.ticks = 3, col.axis = 3)
    mtext("Cross-sectional standard deviation", side = 4, line = 3, col = 3)
}
if(save_plots) dev.off()




## params$D <- 0.38
## X <- sde(xinit, control$times, doublewell, params, control)
## examtimes <- seq(25, 50, by = deltaT)/deltaT
## test <- runner(X[examtimes, -1], sd, k = 1, na_pad = TRUE)

## dev.new(height = 6.5, width = 7)
## par(mai = c(1, 1, 0.5, 1))
## time_ev(X)
## par(new = TRUE)
## plot(examtimes, test, type = "l", col = 2, axes = FALSE, xlab = "", ylab = "")
## axis(4, at = pretty(test), labels = pretty(test), col = 2, col.ticks = 2, col.axis = 2)
## mtext("Cross-sectional standard deviation", 4, 4, col = 2)

## test <- runner(X[, -1], f = function(x)  apply(x, 2, sd), k = 100, na_pad = TRUE)
## test <- do.call(rbind, test[!is.na(test)])

## test <- matrix(rep(1:10, 3), ncol = 3)
## ## testrun <- runner(test, colMeans, 3, na_pad = TRUE)
## testrun <- runner(test, function(x) apply(x, 2, sd), k = 3, na_pad = TRUE)
## testrun <- do.call(rbind, testrun[!is.na(testrun)])
