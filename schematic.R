library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
library(sfsmisc)
source("calc-functions.R")

dynamics <- "doublewell"

g <- make_lattice(length = 10, dim = 2)
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)

model <- get(dynamics)
modelparams <- .doublewell
lmax <- 100
deltaT <- 0.01

bparam <- "u"
xinit <- rep(modelparams$xinit.high, N)

params <- c(modelparams, list(A = A))
control <- list(ncores = ncores, times = 0:50, deltaT = deltaT)

rng <- seq(0, -5, length.out = lmax)

df <- solve_in_range(rng, bparam, model, xinit, params, control, "sde")
attr(df, "model") <- dynamics
attr(df, "network") <- "squarelattice"
attr(df, "simkind") <- "sde"
attr(df, "simtime") <- 0:50
attr(df, "bparam") <- bparam
attr(df, "direction") <- "down"
attr(df, "bparam.vals") <- rng
attr(df, "params") <- params
attr(df, "control") <- control
attr(df, "ntrials") <- 1

moranI <- apply(df, 1, global_moran, A = A)
ssd <- apply(df, 1, sd)

taus <- list(
    moranI = get_tau(df, moranI, TRUE),
    ssd = get_tau(df, ssd, TRUE)
)

farsample <- get_samples(df, "far", n = n)
nearsample <- get_samples(df, "near", n = n)
slopes <- list(
    moranI = list(
        far = get_slope(df, moranI, "far", n = n, return.model = TRUE),
        near = get_slope(df, moranI, "near", n = n, return.model = TRUE)
    ),
    ssd = list(
        far = get_slope(df, ssd, "far", n = n, return.model = TRUE),
        near = get_slope(df, ssd, "near", n = n, return.model = TRUE)
    )
)

tests <- list(
    moranI = classify(df, moranI, n = n),
    ssd = classify(df, ssd, n = n)
)

idx <- get_idx(df)

palette("Okabe-Ito")
dev.new(width = 15, height = 5)
par(mfrow = c(1, 3))
bifplot(df, rng, col = 9, axes = FALSE)
eaxis(1)
eaxis(2)
plot(rng[idx], moranI[idx], xlim = range(rng), axes = FALSE, col = 2, type = "l")
lines(rng[farsample], predict(slopes$moran$far), lwd = 2, lty = 2, col = 1)
lines(rng[nearsample], predict(slopes$moran$near), lwd = 2, lty = 2, col = 1)
eaxis(1)
eaxis(2)
plot(rng[idx], ssd[idx], xlim = range(rng), axes = FALSE, col = 3, type = "l")
lines(rng[farsample], predict(slopes$ssd$far), lwd = 2, lty = 2, col = 1)
lines(rng[nearsample], predict(slopes$ssd$near), lwd = 2, lty = 2, col = 1)
eaxis(1)
eaxis(2)
