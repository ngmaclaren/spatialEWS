load("./data/EWS-data.RData")

library(sfsmisc)
palette("Okabe-Ito")

investigate_lattice <- FALSE

plotdata <- data.frame(
    tau.I = taus$moranI, # these are 'sign-adjusted'
    tau.sd = taus$ssd,
    tau.sk = taus$skew,
    tau.ku = taus$kurt,
    dynamics = dyns,
    network = nets,
    cparam = bparams,
    direction = directions,
    row.names = seq_along(dfs)
)
plotdata$dynamics <- factor(plotdata$dynamics, levels = rev(c("doublewell", "mutualistic", "SIS", "genereg")))
plotdata$network <- factor(plotdata$network)
plotdata$cparam <- factor(plotdata$cparam)
plotdata$direction <- factor(plotdata$direction)

plotdata <- reshape(
    plotdata, varying = c("tau.I", "tau.sd", "tau.sk", "tau.ku"),
    v.names = "tau", timevar = "EWS", times = c("I", "s", "g1", "g2"),
    direction = "long", new.row.names = 1:10000
)

plotdata$EWS <- factor(plotdata$EWS, levels = c("I", "s", "g1", "g2"))
find_adjust <- function(cparam, EWS) {
    cparam <- as.character(cparam)
    EWS <- as.character(EWS)
    if(cparam == "D") {
        ##switch(EWS, I = 0.45, s = 0.32, g1 = 0.19, g2 = 0.06)
        switch(EWS, I = 0.39, s = 0.28, g1 = 0.17, g2 = 0.06)
    } else if(cparam == "u") {
        ##switch(EWS, I = -0.06, s = -0.19, g1 = -0.32, g2 = -0.45)
        switch(EWS, I = -0.06, s = -0.17, g1 = -0.28, g2 = -0.39)
    }
}
plotdata$adjust <- mapply(find_adjust, plotdata$cparam, plotdata$EWS)
plotdata$pch <- ifelse(as.character(plotdata$cparam) == "D", 1, 4)
plotdata$col <- adjustcolor(as.integer(plotdata$EWS) + 1, 1)
plotdata$col <- gsub(palette()[5], palette()[8], plotdata$col)

                                        # For lattice investigation
if(investigate_lattice) {
    classifications <- data.frame(
        I = mapply(classify, dfs, moranIs, n = 5),
        s = mapply(classify, dfs, ssds, n = 5),
        g1 = mapply(classify, dfs, skews, n = 5),
        g2 = mapply(classify, dfs, kurts, n = 5),
        dynamics = dyns,
        network = nets,
        bparam = bparams,
        direction = directions
    )
    rownames(classifications) <- seq_len(nrow(classifications))
    subset(classifications, network == "lattice")

    plotdata$col[plotdata$network != "lattice"] <- gsub("FF", "00", plotdata$col[plotdata$network != "lattice"])
} else {
    plotdata <- subset(plotdata, network != "lattice")
}

ht <- 7; wd <- 12
if(save_plots) {
    if(investigate_lattice) {
        pdf("./img/tauplot-lattice.pdf", height = ht, width = wd)
    } else {
        pdf("./img/tauplot-input.pdf", height = ht, width = wd)
    }
} else {
    dev.new(height = ht, width = wd)
}
labelsize <- 1.75
ptcex <- 1.5
jitter.amount <- 0.01
pt.lwd <- 1.5
par(mfrow = c(1, 2), mar = c(4, 10, 1.5, 0.5))
with(subset(plotdata, direction == "up"), {
    xlim <- c(-1, 1)
    ylim <- c(0.5, 4.5)
    plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    abline(v = 0)
    box()
    axis(2, at = 1:4, labels = levels(dynamics), cex.axis = labelsize, las = 2)
    eaxis(1, cex.axis = labelsize)
    title(xlab = "Sign-adjusted tau", cex.lab = labelsize)
    mtext("(a) Ascending", line = 0.2, adj = 0, cex = labelsize)
    points(
        tau, jitter(as.numeric(dynamics) + adjust, amount = jitter.amount),
        pch = pch, col = col, cex = ptcex, lwd = pt.lwd
    )
    legend("bottomleft", bty = "n", ncol = 1, col = c(2:3, 9, 9), pch = c(1, 1, 1, 4),
           cex = 0.75*labelsize, pt.cex = ptcex, pt.lwd = 2,
           legend = c("Moran's I", "SSD", "D", "u"))
})
with(subset(plotdata, direction == "down"), {
    xlim <- c(-1, 1)
    ylim <- c(0.5, 4.5)
    plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    abline(v = 0)
    box()
    axis(2, at = 1:4, labels = levels(dynamics), cex.axis = labelsize, las = 2)
    eaxis(1, cex.axis = labelsize)
    title(xlab = "Sign-adjusted tau", cex.lab = labelsize)
    mtext("(b) Descending", line = 0.2, adj = 0, cex = labelsize)
    points(
        tau, jitter(as.numeric(dynamics) + adjust, amount = jitter.amount),
        pch = pch, col = col, cex = ptcex, lwd = pt.lwd
    )
})
if(save_plots) dev.off()


                                        # reshape classifications
classified <- reshape(
    classifications,
    varying = c("I", "s", "g1", "g2"),
    v.names = "class",
    timevar = "EWS",
    times = c("I", "s", "g1", "g2"),
    direction = "long"
)
rownames(classified) <- seq(nrow(classified))
colnames(classified)[which(colnames(classified) == "bparam")] <- "cparam"

df <- merge(plotdata, classified, by = c("dynamics", "network", "cparam", "direction", "EWS"))
df <- df[, !grepl("id.", colnames(df))]
rdf <- subset(df, EWS %in% c("I", "s"))
rdf$mainclass <- ifelse(rdf$class == "unsuccessful", "unsuccessful", "successful")

tauagg <- aggregate(tau ~ network + EWS + mainclass, data = rdf, FUN = mean)

subset(tauagg, network == "lattice")
##     network EWS    mainclass       tau
## 20  lattice   I   successful 0.2711190
## 55  lattice   s   successful 0.6876040
## 90  lattice   I unsuccessful 0.2056556
## 125 lattice   s unsuccessful 0.1370968
aggregate(tau ~ EWS + mainclass, data = subset(tauagg, network != "lattice"), FUN = mean)
##   EWS    mainclass         tau
## 1   I   successful  0.85730952
## 2   s   successful  0.86373746
## 3   I unsuccessful -0.04010949
## 4   s unsuccessful -0.65842740

any(subset(rdf, network == "lattice" & mainclass == "successful")$tau < 0) # FALSE
any(subset(rdf, network != "lattice" & mainclass == "successful")$tau < 0) # TRUE

with(list(df = subset(rdf, network != "lattice" & mainclass == "successful")),
     df[df$tau < 0, ])
## all classified as reversing

with(list(df = subset(rdf, network != "lattice" & mainclass == "successful")),
     mean(df[df$tau < 0, "tau"]))
## avg tau = -0.731
