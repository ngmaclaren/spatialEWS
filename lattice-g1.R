library(sfsmisc)
library(igraph)
source("calc-functions.R")
save_plots <- TRUE # FALSE
palette("Okabe-Ito")

filenames <- paste0(
    "./data/sims/",
    grep("lattice", list.files("./data/sims/"), value = TRUE)
)

sims <- lapply(filenames, readRDS)

colors <- list(
    base = 1,
    moran = 2,
    ssd = 3,
    skew = 4,
    kurt = 8,
    xi = 9
)

plotit <- function(sim, ...) {
    df <- promote_df(sim)
    A <- attr(df, "params")$A
    cparam.vals <- attr(df, "bparam.vals")
    cparam.name <- attr(df, "bparam")

    idx <- get_idx(df)

    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew

    lwd <- 3

    par(mar = rep(5, 4))
    matplot(cparam.vals, df, type = "l", lty = 1, lwd = 0.5, col = colors$xi, xlab = "", ylab = "", axes = FALSE)
    eaxis(1, cex.axis = 1.75)
    eaxis(2, cex.axis = 1.75)
    title(xlab = cparam.name, cex.lab = 2)
    title(ylab = expression(x[i]), cex.lab = 2)
    title(main = attr(df, "model"))

    par(new = TRUE)
    plot(cparam.vals[idx], skew[idx], type = "l", xlim = range(cparam.vals), lty = 1, lwd = lwd, col = colors$skew,
         xlab = "", ylab = "", axes = FALSE)
    eaxis(
        4, cex.axis = 1.75, col = colors$skew, col.axis = colors$skew,
        small.args = list(col = colors$skew)
    )
}

get_classification <- function(sim, ...) {
    df <- promote_df(sim)
    skew <- apply(df, 1, moments::skewness)
    if(attr(df, "direction") == "down") skew <- -skew
    classify(df, skew, n = 5)
}


ht <- 10
wd <- 10
if(save_plots) {
    pdf("./img/lattice-g1.pdf", height = ht, width = wd)
    ## png("./img/lattice-g1.png", height = ht, width = wd, units = "in", res = 300)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(4, 3))
for(i in seq_along(sims)) {
    plotit(sims[[i]])
    mtext(paste0("(", LETTERS[i], ")"), cex = 1.5, adj = 0.02)
}
plot(NULL, xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = rev(c(0, 11)), main = "Classification")
text(
    y = seq_along(sims), x = rep(0.5, length(sims)),
    labels = sapply(seq_along(sims), function(i) paste0("(", LETTERS[i], ") ", get_classification(sims[[i]])))
)
if(save_plots) dev.off()
