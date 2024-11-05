library(sfsmisc)
library(igraph)
source("calc-functions.R")
save_plots <- FALSE # TRUE
palette("Okabe-Ito")

filenames <- paste0(
    "./data/sims/",
    grep("lattice", list.files("./data/sims/"), value = TRUE)
)
i.dw <- c(4, 3, 2, 1) # dw: u up, u down, D up, D down,
i.ms <- c(8, 7)       # ms: u down, D down,
i.sis <- c(10, 9)     # sis: D up, D down
i.gr <- c(6, 5)       # gr: u down, D down
filenames <- filenames[c(i.dw, i.ms, i.sis, i.gr)]

sims <- lapply(filenames, readRDS)
sims <- lapply(sims, promote_df)

ylims.xi <- list(
    doublewell = range(unlist(sims[1:4])),
    mutualistic = range(unlist(sims[5:6])),
    SIS = range(unlist(sims[7:8])),
    genereg = range(unlist(sims[9:10]))
)

skews <- lapply(sims, function(df) {
    idx <- get_idx(df)
    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew
    skew[idx]
})
ylims.g1 <- list(
    doublewell = range(unlist(skews[1:4])),
    mutualistic = range(unlist(skews[5:6])),
    SIS = range(unlist(skews[7:8])),
    genereg = range(unlist(skews[9:10]))
)

colors <- list(
    base = 1,
    moran = 2,
    ssd = 3,
    skew = 4,
    kurt = 8,
    xi = 9
)

plotit <- function(df, ...) {
    ## df <- promote_df(sim)
    A <- attr(df, "params")$A
    cparam.vals <- attr(df, "bparam.vals")
    cparam.name <- attr(df, "bparam")
    model <- attr(df, "model")

    idx <- get_idx(df)

    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew

    ylim.xi <- ylims.xi[[model]]
    ylim.g1 <- ylims.g1[[model]]
    lwd <- 3

    par(mar = rep(5, 4))
    matplot(cparam.vals, df, type = "l", lty = 1, lwd = 0.5, col = colors$xi, ylim = ylim.xi,
            xlab = "", ylab = "", axes = FALSE)
    eaxis(1, cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    title(xlab = cparam.name, cex.lab = 2)
    title(ylab = expression(x[i]), cex.lab = 2)
    title(main = model)

    par(new = TRUE)
    plot(cparam.vals[idx], skew[idx], type = "l", xlim = range(cparam.vals), lty = 1, lwd = lwd, col = colors$skew,
         ylim = ylim.g1, xlab = "", ylab = "", axes = FALSE)
    eaxis(
        4, cex.axis = 1.75, col = colors$skew, col.axis = colors$skew,
        small.args = list(col = colors$skew)
    )
}

get_classification <- function(df, ...) {
    ## df <- promote_df(sim)
    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew
    classify(df, skew, n = 5)
}
classifications <- sapply(sims, get_classification)

wd <- 3.5*5
ht <- 3.5*2
if(save_plots) {
    pdf("./img/lattice-g1.pdf", height = ht, width = wd)
    ## png("./img/lattice-g1.png", height = ht, width = wd, units = "in", res = 300)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(2, 5))
for(i in seq_along(sims)) {
    plotit(sims[[i]])
    mtext(paste0("(", letters[i], ")"), cex = 1.25, adj = 0.02)
    mtext(tools::toTitleCase(paste(classifications[i])), col = colors$skew, cex = 1.25, adj = 0.02, line = -2)
}
## plot(NULL, xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = rev(c(0, 11)), main = "Classification")
## text(
##     y = seq_along(sims), x = rep(0.5, length(sims)),
##     labels = sapply(seq_along(sims), function(i) paste0("(", LETTERS[i], ") ", get_classification(sims[[i]])))
## )
if(save_plots) dev.off()
