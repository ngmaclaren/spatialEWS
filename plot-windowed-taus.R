
library(parallel)
ncores <- detectCores()-1
source("calc-functions.R")
source("plot-functions.R")
save_plots <- !interactive()

check <- function(df) cat(show(df[1:5, 1:5]), show(df[(nrow(df)-4):nrow(df), 1:5]))

networks <- readRDS("networks.rds")

basin <- list(
    doublewell = 3,
    genereg = 0.1,
    mutualistic = 1,
    SIS = 0.1
)
reduce <- function(dl) { # for now, don't need the extra simulations, just take the first one
    attrs <- attributes(dl)
    df <- dl[[1]]
    attr(df, "info") <- attrs
    return(df)
}
focus <- function(df) { # reduce each df to just the data before the first node transition
    info <- attr(df, "info")
    idx <- switch(
        info$direction,
        down = which(apply(df, 1, min) > basin[[info$model]]),
        up = which(apply(df, 1, max) < basin[[info$model]])
    )
    if(info$direction == "down" & info$bparam == "D") idx <- rev(idx)
    newdf <- df[idx, ]
    info$bparam.vals <- info$bparam.vals[idx]
    attr(newdf, "info") <- info
    return(newdf)
}
summarize <- function(df) { # compute the raw data we need: both metrics and the mean state
    info <- attr(df, "info")
    moranI <- apply(df, 1, global_moran, info$params$A)
    ssd <- apply(df, 1, sd)
    meanstate <- rowMeans(df)
    newdf <- data.frame(moranI = moranI, ssd = ssd, meanstate = meanstate, bparam = info$bparam.vals)
    attr(newdf, "info") <- info
    return(newdf)
}

datadir <- "./data/"
datafiles <- list.files(datadir, pattern = "*.rds")

alldata <- mcMap(readRDS, paste0(datadir, datafiles), mc.cores = ncores)
## need to fix the names, for convenience
dfs <- lapply(alldata, reduce)
qvars <- list(
    model = sapply(dfs, function(df) attr(df, "info")$model),
    network = sapply(dfs, function(df) attr(df, "info")$network),
    direction = sapply(dfs, function(df) attr(df, "info")$direction),
    bparam = sapply(dfs, function(df) attr(df, "info")$bparam)
)
                                        # for debugging
## thisone <- which(
##     qvars$model == "doublewell" & qvars$network == "ba" & qvars$direction == "down" & qvars$bparam == "D"
## )
## check(dfs[[thisone]])
## When bparam is D and direction is down, need to reverse the order of the bparam.vals so that we go from large to small. Can fix that in the simulation eventually, but for this morning just reverse it manually here.
dfs <- mclapply(dfs, focus, mc.cores = ncores)
dat <- mclapply(dfs, summarize, mc.cores = ncores) # is the order always correct?
rm(alldata)

## ??
## table(qvars$model)
## doublewell     genereg mutualistic         SIS 
##        167         168         167          83 
## table(qvars$network)

##    akatoreA    akatoreB          ba     berwick   blackrock     catlins  chesapeake    coweeta1   coweeta17 
##          14          14          14          14          14          14          14          14          14 
##  dempstersA  dempstersS dempstersSu     dolphin          er          ff  genefusion      german         ggs 
##          14          14          14          14          14          14          14          14          14 
##         ggt       healy        jazz     kyeburn  lilkyeburn     martins  narrowdale      netsci    northcol 
##          14          14          14          14          14          14          14          11          14 
##     pdzbase          pl          sl       stony    suttonAu    suttonSp    suttonSu         sw3         sw4 
##          14          14          14          14          14          14          14          14          14 
##         sw5          t3          t4          t5        troy      venlaw 
##          14          14          14          14          14          14
## table(qvars$direction)

## down   up 
##  293  292 

## compute the sequence of tau values (at fifth to last, 10th to last, etc.)
window_taus <- function(df, metric, window.length = 5) { # better do a match.arg for metric
    info <- attr(df, "info")
    end <- nrow(df)
    if(end < window.length) return(NA)
    start <- end - (window.length-1)
    windows <- rev(seq(start, 1, by = -window.length)) # the start of each window, from longest to shortest (i.e., closest to the first node transition
    taus <- sapply(
        windows,
        function(w) {
            x <- df[, metric][w:end]
            y <- df[, "bparam"][w:end]
            cor(x, y, method = "kendall")
        }
    )
    windowlength <- (end - windows) + 1
    newdf <- data.frame(tau = taus, windowlength = windowlength)
    return(newdf)
}

                                        # 49 NA values for each of the below. Why?
taus.moranI <- mclapply(dat, FUN = window_taus, "moranI", mc.cores = ncores)
taus.ssd <- mclapply(dat, FUN = window_taus, "ssd", mc.cores = ncores)

## make the plot frame(s)
selected <- which(# across networks, fix model, direction, and bparam for now. Might want to assign features to some
    qvars$model == "mutualistic" &
    qvars$direction == "down" &
    qvars$bparam == "D"
)
plotdata <- list(
    moranI = taus.moranI[selected],
    ssd = taus.ssd[selected]
)
plotinfo <- lapply(qvars, function(x) x[selected])
xlim <- range(
    c(unlist(lapply(plotdata$moranI, function(df) df$windowlength)),
      unlist(lapply(plotdata$ssd, function(df) df$windowlength))
      )
)
ylim <- c(-1, 1)
palette("Set 1")
plot(NULL, xlim = rev(xlim), ylim = ylim, xlab = "Window length", ylab = expression(tau))
abline(h = 0, lwd = .5, lty = 3)
## make the lines. Should have e.g. 167/2(?) doublewell-up lines.
for(i in seq_along(selected)) {
    lines(plotdata$moranI[[i]]$windowlength, plotdata$moranI[[i]]$tau, lty = 1, lwd = 0.5, col = 1)
    lines(plotdata$ssd[[i]]$windowlength, plotdata$ssd[[i]]$tau, lty = 1, lwd = 0.5, col = 2)
}
legend(
    x = xlim[2], y = 0.5, bty = "n",
    col = 1:2, lwd = 2,
    legend = c("Moran's I", "Std. Dev")
)
