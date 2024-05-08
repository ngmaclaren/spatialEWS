## Working on classify().
## Is it basically working? I think so, but need to figure out to check further.
## Certainly the parameters, multiple and theta, must be bad. Musn't they?

library(igraph)
library(localsolver)
source("calc-functions.R")
palette("Okabe-Ito")

networks <- readRDS("./data/networks.rds")

simdir <- "./data/sims/"
simfiles <- list.files(simdir)
simdata <- lapply(paste0(simdir, simfiles), function(simfile) readRDS(simfile))

dyns <- sapply(simdata, attr, "model")
nets <- sapply(simdata, attr, "network")
bparams <- sapply(simdata, attr, "bparam")
directions <- sapply(simdata, attr, "direction")
cparam.vals <- lapply(simdata, attr, "bparam.vals")
As <- lapply(simdata, function(dat) attr(dat, "params")$A)

## Continue here instead ##
dfs <- lapply(
    simdata,
    function(sim) {
        attribs <- attributes(sim)
        df <- sim[[1]]
                                        #
        ## if(attribs$direction == "down") {
        ##     df <- df[rev(seq(nrow(df))), ]
        ##     attribs$bparam.vals <- rev(attribs$bparam.vals)
        ## }
                                        #
        attributes(df) <- c(attributes(df), attribs)
        df
    }
)


dfnames <- mapply(
    function(dyn, net, bp, dir) paste(dyn, net, bp, dir, sep = "-"),
    dyns, nets, bparams, directions,
    USE.NAMES = FALSE
)
names(dfs) <- dfnames

u.down <- which(grepl("-u-", dfnames, fixed = TRUE) & grepl("-down", dfnames, fixed = TRUE))
u.up <- which(grepl("-u-", dfnames, fixed = TRUE) & grepl("-up", dfnames, fixed = TRUE))
D.down <- which(grepl("-D-", dfnames, fixed = TRUE) & grepl("-down", dfnames, fixed = TRUE))
D.up <- which(grepl("-D-", dfnames, fixed = TRUE) & grepl("-up", dfnames, fixed = TRUE))

moranIs<- mapply(
    function(df, A) apply(df, 1, global_moran, A = A),
    dfs, As, SIMPLIFY = FALSE, USE.NAMES = TRUE
)
cssds <- lapply(dfs, function(df) apply(df, 1, sd))

basins <- list( ## GLOBAL
    doublewell = 3,
    genereg = 0.1,
    mutualistic = 1,
    SIS = 0.1
)

## Now, all the information I need should be with the df
get_idx <- function(df) {
    direction <- attr(df, "direction")
    dynamics <- attr(df, "model")
    basin <- basins[[dynamics]]

    switch(
        direction,
        down = which(apply(df, 1, min) > basin),
        up = which(apply(df, 1, max) < basin)
    )
}

get_samples <- function(df, which = c("near", "far"), n = 5) {
    direction <- attr(df, "direction")
    whichslope <- match.arg(which)

    idx <- get_idx(df)

    switch(
        direction,
        up = switch(whichslope, far = idx[1:n], near = rev(idx[seq(length(idx), by = -1, length.out = n)])),
        down = switch(whichslope, near = idx[1:n], far = rev(idx[seq(length(idx), by = -1, length.out = n)]))
    )
}

## get_quadratic_fit <- function(EWS, cparam.vals, direction) {}

get_first_slope <- function(df, EWS, return.model = FALSE) { # compute the EWS ahead of time? Yes.
    direction <- attr(df, "direction")
    
    samples <- get_samples(df, "far") # expose n?

    cparam <- attr(df, "bparam.vals")[samples]
    ews <- EWS[samples]

    m <- lm(ews ~ cparam)
    
    if(return.model) {
        return(m)
    } else {
        return(as.numeric(coef(m)[2]))
    }
}
    
get_second_slope <- function(df, EWS, return.model = FALSE) { # compute the EWS ahead of time? Yes.
    direction <- attr(df, "direction")
    
    samples <- get_samples(df, "near") # expose n?

    cparam <- attr(df, "bparam.vals")[samples]
    ews <- EWS[samples]

    m <- lm(ews ~ cparam)

    if(return.model) {
        return(m)
    } else {
        return(as.numeric(coef(m)[2]))
    }
}

get_quadratic_fit <- function(df, EWS) {
                                        # I think this function will be no good: it imposes too many assumptions
                                        # on the data.
                                        # If used, negative quad term means arms of parabola point down
    idx <- get_idx(df)

    cparam <- attr(df, "bparam.vals")[idx]
    ews <- EWS[idx]

    fit <- lm(ews ~ poly(cparam, 2))
    return(fit)
}

dev.new(height = 10, width = 10)
par(mfrow = c(2, 2))
plot(
    mapply(get_first_slope, dfs, cssds)[D.up], mapply(get_second_slope, dfs, cssds)[D.up],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "D up"
)
mtext("Cross-sectional standard deviation", adj = 1, line = 3, cex = 1, font = 2)
plot(
    mapply(get_first_slope, dfs, cssds)[D.down], mapply(get_second_slope, dfs, cssds)[D.down],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "D down"
)
plot(
    mapply(get_first_slope, dfs, cssds)[u.up], mapply(get_second_slope, dfs, cssds)[u.up],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "u up"
)
plot(
    mapply(get_first_slope, dfs, cssds)[u.down], mapply(get_second_slope, dfs, cssds)[u.down],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "u down"
)

dev.new(height = 10, width = 10)
par(mfrow = c(2, 2))
plot(
    mapply(get_first_slope, dfs, moranIs)[D.up], mapply(get_second_slope, dfs, moranIs)[D.up],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "D up"
)
mtext("Moran's I", adj = 1, line = 3, cex = 1, font = 2)
plot(
    mapply(get_first_slope, dfs, moranIs)[D.down], mapply(get_second_slope, dfs, moranIs)[D.down],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "D down"
)
plot(
    mapply(get_first_slope, dfs, moranIs)[u.up], mapply(get_second_slope, dfs, moranIs)[u.up],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "u up"
)
plot(
    mapply(get_first_slope, dfs, moranIs)[u.down], mapply(get_second_slope, dfs, moranIs)[u.down],
    xlab = "Far from the bifurcation", ylab = "Near the bifurcation",
    cex = 1.5, pch = 16, cex.axis = 1.5, cex.lab = 1.5, main = "u down"
)

f.slope <- data.frame(
    cssd = mapply(get_first_slope, dfs, cssds),
    moranI = mapply(get_first_slope, dfs, moranIs),
    dynamics = dyns,
    network = nets,
    bparam = bparams,
    direction = directions,
    slope = "first"
)
s.slope <- data.frame(
    cssd = mapply(get_second_slope, dfs, cssds),
    moranI = mapply(get_second_slope, dfs, moranIs),
    dynamics = dyns,
    network = nets,
    bparam = bparams,
    direction = directions,
    slope = "second"
)
slopes <- rbind(f.slope, s.slope)
hist(slopes$cssd[slopes$slope == "second"])


classify <- function(df, EWS, multiple = 5, threshold = 5) {
    firstslope <- get_first_slope(df, EWS)
    secondslope <- get_second_slope(df, EWS)

    firstsign <- sign(firstslope)
    secondsign <- sign(secondslope)
    mult <- secondslope/firstslope

    direction <- attr(df, "direction")

    if(direction == "up") {
        if(sign(firstslope) > 0 & sign(secondslope) > 0 & abs(mult) > multiple & secondslope > threshold) {
            return("consistent")
        } else if(sign(firstslope) < 0 & sign(secondslope) > 0 & abs(mult) > multiple & secondslope > threshold) {
            return("signswitch")
        } else {
            return("bad")
        }
    } else { # down
        if(sign(firstslope) < 0 & sign(secondslope) < 0 & abs(mult) > multiple & secondslope > threshold) {
            return("consistent")
        } else if(sign(firstslope) > 0 & sign(secondslope) < 0 & abs(mult) > multiple & secondslope > threshold) {
            return("signswitch")
        } else {
            return("bad")
        }
    }
}


## checking...
for(i in 1:4) {
    df <- dfs[[i]]
    EWS <- cssds[[i]]
    idx <- get_idx(df)
    cparam <- attr(df, "bparam.vals")[idx]
    ews <- EWS[idx]
    nearsamples <- get_samples(df, "near")
    farsamples <- get_samples(df, "far")
    fsm <- get_first_slope(df, EWS, TRUE)
    ssm <- get_second_slope(df, EWS, TRUE)
    qf <- get_quadratic_fit(df, EWS)
    dev.new()
    palette("R4")
    plot(cparam, ews, xlab = attr(df, "bparam"), ylab = "Moran's I")
    lines(cparam, predict(qf), lty = 2, lwd = 1, col = 2)
    lines(attr(df, "bparam.vals")[nearsamples], predict(ssm), lty = 1, lwd = 2, col = 3)
    lines(attr(df, "bparam.vals")[farsamples], predict(fsm), lty = 1, lwd = 2, col = 3)
    palette("Okabe-Ito")
}




pdf("./img/diagnostic-plots.pdf")
for(i in seq_along(dfs)) {
    if(isFALSE(i %in% u.down)) next
    
    bifplot(dfs[[i]], cparam.vals[[i]], col = adjustcolor(1, 0.5), lwd = 0.5, main = dfnames[i])

    basin <- basins[[dyns[i]]]
    idx <- switch(
        directions[i],
        down = which(apply(dfs[[i]], 1, min) > basin),
        up = which(apply(dfs[[i]], 1, max) < basin)
    )

    par(new = TRUE)
    plot(
        cparam.vals[[i]][idx], moranIs[[i]][idx], type = "l", col = 2, lwd = 2, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals[[i]])
    )
    mtext(paste("Moran's I: first =", round(f.slope$moranI[i], 3),
                ", second =", round(s.slope$moranI[i], 3)),
          line = -1)
    mtext(paste("CSSD: first =", round(f.slope$cssd[i], 3),
                ", second =", round(s.slope$cssd[i], 3)),
          line = -3)
    par(new = TRUE)
    plot(
        cparam.vals[[i]][idx], cssds[[i]][idx], type = "l", col = 3, lwd = 2, lty = 1,
        axes = FALSE, xlab = "", ylab = "", xlim = range(cparam.vals[[i]])
    )

    legend("bottomright", bty = "n", lwd = 2, col = 2:3, lty = 1,
           legend = c("Moran's I", "Cross-sectional SD"))
    
}
dev.off()
