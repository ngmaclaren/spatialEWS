                                        # Generates the raw data for Fig 6

### Use the stored BA network sims
### Use u as the bifurcation parameter
### Show down, with arrows, and up, with arrows.
### Highlight some things.
library(sfsmisc)
source("calc-functions.R")

sims <- list(
    down = readRDS("./data/sims/doublewell-barabasialbert-u-sde-50-down-0-1-TRUE.rds"),
    up = readRDS("./data/sims/doublewell-barabasialbert-u-sde-50-up-0-1-TRUE.rds")
)
sims <- lapply(sims, promote_df)

palettecolors <- c(blue = "#3584e4", green = "#33d17a", yellow = "#f6d32d", orange = "#ff7800", red = "#e01b24",
                   purple = "#9141ac", brown = "#986a44", light = "#deddda", dark = "#3d3846")
palette(palettecolors)
labelsize <- 1.5

pdf("./img/schematic-panels.pdf", width = 21, height = 7)
par(mfrow = c(1, 3))
with(list(df = sims$down), {
    par(mai = c(1.25, 1.25, 0.75, 0.75))
    matplot(
        attr(df, "bparam.vals"), df, type = "l", lty = 1, lwd = 0.25, col = "black",
        axes = FALSE, xlab = "u", ylab = expression(x[i]), cex.lab = labelsize
    )
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
})

with(list(df = sims$up), {
    par(mai = c(1.25, 1.25, 0.75, 0.75))
    matplot(
        attr(df, "bparam.vals"), df, type = "l", lty = 1, lwd = 0.25, col = "black",
        axes = FALSE, xlab = "u", ylab = expression(x[i]), cex.lab = labelsize
    )
    eaxis(1, cex.axis = labelsize)
    eaxis(2, cex.axis = labelsize)
})

with(list(df = sims$up), {
    bparam <- attr(df, "bparam.vals")
    moranI <- apply(df, 1, global_moran, A = attr(df, "params")$A)
    ssd <- apply(df, 1, sd)
    idx <- get_idx(df)
    par(mai = c(1.25, 1.25, 0.75, 0.75))
    matplot(
        bparam, df, type = "l", lty = 1, lwd = 0.25, col = 8,
        axes = FALSE, xlab = "", ylab = "", cex.lab = labelsize,
        ylim = c(min(df), 3)
    )
    par(new = TRUE)
    plot(
        bparam[idx], moranI[idx], type = "l", lty = 1, lwd = 2, col = palettecolors["blue"], xlim = range(bparam),
        axes = FALSE, xlab = "", ylab = "", cex.lab = labelsize
    )
    eaxis(1, cex.axis = labelsize)
    title(xlab = "u", cex.lab = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(ylab = "Early warning signal", cex.lab = labelsize, line = 4)
})
dev.off()


pdf("./img/anotherpanel.pdf")
with(list(df = sims$down), {
    bparam <- attr(df, "bparam.vals")
    moranI <- apply(df, 1, global_moran, A = attr(df, "params")$A)
    ssd <- apply(df, 1, sd)
    idx <- get_idx(df)
    par(mai = c(1.25, 1.25, 0.75, 0.75))
    matplot(
        bparam, df, type = "l", lty = 1, lwd = 0.25, col = 8,
        axes = FALSE, xlab = "", ylab = "", cex.lab = labelsize,
        ylim = c(4, max(df))
    )
    par(new = TRUE)
    plot(
        bparam[idx], moranI[idx], type = "l", lty = 1, lwd = 2, col = palettecolors["blue"], xlim = range(bparam),
        axes = FALSE, xlab = "", ylab = ""
    )
    eaxis(1, cex.axis = labelsize)
    title(xlab = "u", cex.lab = labelsize)
    eaxis(2, cex.axis = labelsize)
    title(ylab = "Early warning signal", cex.lab = labelsize, line = 4)
})
dev.off()

df <- promote_df(readRDS("./data/sims/doublewell-barabasialbert-u-sde-50-down-0-1-TRUE.rds"))
moranI <- apply(df, 1, global_moran, A = attr(df, "params")$A)
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
