global_moran <- function(x, A) {# rewrite to accept adjacency matrix instead of network
    N <- length(x)
    stopifnot(N == nrow(A))
    stopifnot(nrow(A) == ncol(A))

    W <- sum(A)
    x <- x - mean(x) # center

    (N/W)*sum(A*outer(x, x))/sum(x^2)
}

local_moran <- function(i, x, A) {
                                        # this is too large by exactly the mean degree
                                        # Wikipedia says I = sum(I/N) but I get I = mean(k)*sum(I/N)
    N <- length(x)
    stopifnot(N == nrow(A))
    stopifnot(nrow(A) == ncol(A))
    x <- as.numeric(x)
    x <- x - mean(x)

    m2 <- sum(x^2)/N
    (x[i]/m2)*sum(A[i, ]*x)
}

## check_zerocrossing <- function() {}
## q1 <- function() {}

## check_magnitude <- function() {}
## q2 <- function() {}

### Classify the EWS
basins <- list( ## GLOBAL
    doublewell = 3, # separatrix
    genereg = 0.005, # 5*(noise strength)
    mutualistic = 1, # Allee constant
    SIS = 0.005 # 5*(noise strength)
)

promote_df <- function(sim, which = 1) {
    attribs <- attributes(sim)
    df <- sim[[which]]
    attributes(df) <- c(attributes(df), attribs)
    df
}

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

get_samples <- function(df, which = c("near", "far"), n = 3) {
    whichslope <- match.arg(which)

    idx <- get_idx(df)

    switch(
        whichslope,
        near = rev(seq(length(idx), by = -1, length.out = n)),
        far = seq(1, by = 1, length.out = n)
    )
}

get_tau <- function(df, EWS, adjust.sign = FALSE) {
    idx <- get_idx(df)
    midpoint <- floor(median(idx))
    idx <- idx[which(idx > midpoint)]

    cparam <- attr(df, "bparam.vals")[idx]
    ews <- EWS[idx]

    tau <- cor(cparam, ews, method = "kendall")

    if(adjust.sign) {
        if(attr(df, "direction") == "down") {
            tau <- -tau
        }
    }

    return(tau)
}

get_slope <- function(df, EWS, which = c("near", "far"), n = 3, return.model = FALSE) {
    whichslope <- match.arg(which, c("near", "far"))

    samples <- get_samples(df, whichslope, n)

    cparam <- attr(df, "bparam.vals")[samples]
    ews <- EWS[samples]
    m <- lm(ews ~ cparam)

    if(return.model) {
        return(m)
    } else {
        return(as.numeric(coef(m)[2]))
    }
}

classify <- function(df, EWS, threshold = 2, n = 3) {
    firstslope <- get_slope(df, EWS, "far", n = n)
    secondslope <- get_slope(df, EWS, "near", n = n)

    firstsign <- sign(firstslope)
    secondsign <- sign(secondslope)
    ratio <- secondslope/firstslope

    direction <- attr(df, "direction")

    if(direction == "up") {
        if(sign(firstslope) > 0 & sign(secondslope) > 0 & abs(ratio) > threshold) {
            return("consistent")
        } else if(sign(firstslope) < 0 & sign(secondslope) > 0 & abs(ratio) > 0.5*threshold) {
            return("signswitch")
        } else {
            return("bad")
        }
    } else { # down
        if(sign(firstslope) < 0 & sign(secondslope) < 0 & abs(ratio) > threshold) {
            return("consistent")
        } else if(sign(firstslope) > 0 & sign(secondslope) < 0 & abs(ratio) > 0.5*threshold) {
            return("signswitch")
        } else {
            return("bad")
        }
    }
}
