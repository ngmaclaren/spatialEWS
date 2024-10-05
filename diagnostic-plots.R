library(sdn)
library(sfsmisc)

load("./data/EWS-data.RData")

diagnostic_plot <- function(df, moranI, ssd, skew, kurt, n = 3, main = "", ...) {
    cparam.vals <- attr(df, "bparam.vals")
    idx <- get_idx(df)
    farsample <- get_samples(df, "far", n = n)
    nearsample <- get_samples(df, "near", n = n)

    farmodel.moran <- get_slope(df, moranI, "far", n = n, return.model = TRUE)
    farmodel.ssd <- get_slope(df, ssd, "far", n = n, return.model = TRUE)
    nearmodel.moran <- get_slope(df, moranI, "near", n = n, return.model = TRUE)
    nearmodel.ssd <- get_slope(df, ssd, "near", n = n, return.model = TRUE)

    test.moran <- classify(df, moranI, n = n)
    test.ssd <- classify(df, ssd, n = n)

    bifplot(df, cparam.vals, col = adjustcolor(1, 0.5), lwd = 0.5, ...) #  ylim = c(0, 0.02),

    mtext("Test results", line = 2, adj = 0, font = 2)
    mtext(paste("Moran's I:", test.moran), line = 1, adj = 0)
    mtext(paste("Spatial SD:", test.ssd), line = 0, adj = 0)

    mtext("Kendall's tau", line = 2, adj = 0.5, font = 2)
    mtext(paste("Moran's I:", round(get_tau(df, moranI, TRUE), 2)), line = 1, adj = 0.5)
    mtext(paste("Spatial SD:", round(get_tau(df, ssd, TRUE), 2)), line = 0, adj = 0.5)

    mtext(paste("Dynamics:", attr(df, "model")), line = 3, adj = 1)
    mtext(paste("Control parameter:", attr(df, "bparam")), line = 2, adj = 1)
    mtext(paste("Direction:", attr(df, "direction")), line = 1, adj = 1)
    mtext(paste("Network:", attr(df, "network")), line = 0, adj = 1)

    par(new = TRUE)
    plot(
        cparam.vals[idx], moranI[idx], type = "l", col = 2, lwd = 3, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    lines(cparam.vals[farsample], predict(farmodel.moran), lwd = 1.5, lty = 1, col = 1)
    lines(cparam.vals[nearsample], predict(nearmodel.moran), lwd = 1.5, lty = 1, col = 1)

    par(new = TRUE)
    plot(
        cparam.vals[idx], ssd[idx], type = "l", col = 3, lwd = 3, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    lines(cparam.vals[farsample], predict(farmodel.ssd), lwd = 1.5, lty = 1, col = 1)
    lines(cparam.vals[nearsample], predict(nearmodel.ssd), lwd = 1.5, lty = 1, col = 1)

    par(new = TRUE)
    plot(
        cparam.vals[idx], skew[idx], type = "l", col = 4, lwd = 3, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )
    eaxis(4, col = 4, col.axis = 4, small.args = list(col = 4))
    par(new = TRUE)
    plot(
        cparam.vals[idx], kurt[idx], type = "l", col = 5, lwd = 3, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals)
    )

    legend("bottomright", bty = "n", lwd = 2, col = 2:5, lty = 1,
           legend = c("Moran's I", "Spatial std. dev.", "Skewness", "Kurtosis"))
}

pdf("./img/diagnostic-plots.pdf")
par(mar = rep(4, 4))
mapply(
    diagnostic_plot, dfs[plotorder], moranIs[plotorder], ssds[plotorder], skews[plotorder], kurts[plotorder], n = n
)
## with(
##     list(theseplots = grep("genereg", dfnames)), 
##     mapply(diagnostic_plot, dfs[theseplots], moranIs[theseplots], ssds[theseplots], n = n)
## )
dev.off()
