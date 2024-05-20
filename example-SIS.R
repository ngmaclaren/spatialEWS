library(igraph)
source("calc-functions.R")

save_plots <- FALSE # TRUE
palette("Okabe-Ito")

g <- readRDS("./data/networks.rds")[["drug"]]
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)

SIS.Du <- readRDS("./data/sims/SIS-drug-D-sde-50-up-0-1-TRUE.rds")
SIS.Dd <- readRDS("./data/sims/SIS-drug-D-sde-50-down-0-1-TRUE.rds")

ht <- 15/2
wd <- 15
if(save_plots) {
    pdf("./img/example-SIS.pdf", height = ht, width = wd)
} else { 
    dev.new(height = ht, width = wd)
}
par(mfrow = c(1, 2), mai = c(1, 1, 0.1, 0.1))
with(
    list(
        X = SIS.Du[[1]],
        D = attr(SIS.Du, "bparam.vals")
    ), {
        maxs <- apply(X, 1, max)
        far1 <- 10
        far2 <- 35
        nearest <- which(maxs < 3)
        nearest <- nearest[length(nearest)]
        near <- nearest - 3
                                        # marking lines; will need the SSDs for these lines
        markD <- D[c(far1, far2, near, nearest)]
        markX <- X[c(far1, far2, near, nearest), ]
        print(apply(markX, 1, sd))
        print(apply(markX, 1, global_moran, A))
                                        # plot
        matplot(
            D, X, type = "l", lty = 1, lwd = 0.5, col = 9,
            ylim = c(0, 0.01), xlim = c(0, 0.2),
            xlab = "D", ylab = expression(x[i]), axes = FALSE, font.lab = 3, cex.lab = 2
        )
        axis(1, cex.axis = 2)
        axis(2, cex.axis = 2)
        segments(
            x0 = markD,
            y0 = apply(markX, 1, min) - 0.005,
            y1 = apply(markX, 1, max) + 0.005,
            col = 3, lty = 1, lwd = 4
        )
    }
)
