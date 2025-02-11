global_moran <- function(x, A) {
    N <- length(x)
    stopifnot(N == nrow(A))
    stopifnot(nrow(A) == ncol(A))

    W <- sum(A)
    x <- x - mean(x) # center, for convenience below

    (N/W)*sum(A*outer(x, x))/sum(x^2)
}

CV <- function(x) {
    sd(x)/mean(x)
}

### Classify the EWS
basins <- list( ## GLOBAL
    doublewell = 3, # separatrix
    genereg = 0.005, # 5*(noise strength), 10*sigma is too large
    mutualistic = 1, # Allee constant
    SIS = 0.005 # 5*(noise strength), same as for genereg
)

                                        # For convenience.
                                        # The simulations can produce lists of replicated experiments, but we only
                                        # use the first (and only) data frame in each list. This function attaches
                                        # all the simulation information, stored as attributes, to the data frame.
promote_df <- function(sim, which = 1) {
    attribs <- attributes(sim)
    df <- sim[[which]]
    attributes(df) <- c(attributes(df), attribs)
    df
}

                                        # This function returns the row indices of the simulation range.
                                        # It implements a decision rule about states which says, "I'm in my original
                                        # state unless I am as or more extreme than a certain value."
get_idx <- function(df, restrict = NULL) {
    direction <- attr(df, "direction")
    dynamics <- attr(df, "model")
    basin <- basins[[dynamics]]

    idx <- switch(
        direction,
        down = which(apply(df, 1, min) > basin),
        up = which(apply(df, 1, max) < basin) 
    )

    if(is.null(restrict)) {
        return(idx)
    } else {
                                        # pass a list of quantiles to restrict
        stopifnot("far" %in% names(restrict) | "near" %in% names(restrict))
        if("far" %in% names(restrict)) far <- restrict$far else far <- 0
        if("near" %in% names(restrict)) near <- restrict$near else near <- 1
        far <- floor(quantile(idx, probs = far))
        near <- ceiling(quantile(idx, probs = near))
        ## idx <- idx[far:length(idx)]
        idx <- idx[far:near]
        return(idx)
    }
}

                                        # Returns the row indices of the first n (far) or last n (near) samples in
                                        # the relevant range.
get_samples <- function(df, which = c("near", "far"), n = 5, restrict = NULL) {
    whichslope <- match.arg(which)

    idx <- get_idx(df, restrict = restrict)

    switch(
        whichslope,
        ## near = rev(seq(length(idx), by = -1, length.out = n)),
        near = rev(idx[seq(length(idx), by = -1, length.out = n)]),
        ## far = seq(1, by = 1, length.out = n)
        far = idx[1:n]
    )
}

                                        # returns the Kendall's tau between the control parameter and the EWS.
                                        # The EWS should be computed beforehand. If up, then τ, else τ' (= -τ).
get_tau <- function(df, EWS, adjust.sign = FALSE, restrict = NULL) {
    idx <- get_idx(df, restrict = restrict)
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

                                        # Returns the simple slope for the first/last n data in the relevant range.
                                        # Independent: control parameter. Dependent: EWS. 
                                        # Optionally returns the whole model.
get_slope <- function(df, EWS, which = c("near", "far"), n = 5, return.model = FALSE, restrict = NULL) {
    whichslope <- match.arg(which, c("near", "far"))

    samples <- get_samples(df, which = whichslope, n = n, restrict = restrict)

    cparam <- attr(df, "bparam.vals")[samples]
    ews <- EWS[samples]
    m <- lm(ews ~ cparam)

    if(return.model) {
        return(m)
    } else {
        return(as.numeric(coef(m)[2]))
    }
}

                                        # Classifies an EWS as accelerating, reversing, or unsuccessful based on
                                        # the slope of the EWS vs. the control parameter and a threshold value
                                        # (default = 2). 
classify <- function(df, EWS, threshold = 2, n = 5, restrict = NULL) {
    firstslope <- get_slope(df, EWS, "far", n = n, restrict = restrict)
    secondslope <- get_slope(df, EWS, "near", n = n, restrict = restrict)

    firstsign <- sign(firstslope)
    secondsign <- sign(secondslope)
    ratio <- secondslope/firstslope

    direction <- attr(df, "direction")

    if(direction == "up") {
        if(sign(firstslope) > 0 & sign(secondslope) > 0 & abs(ratio) > threshold) {
            return("accelerating")
        } else if(sign(firstslope) < 0 & sign(secondslope) > 0 & abs(ratio) > 0.5*threshold) {
            return("reversing")
        } else {
            return("unsuccessful")
        }
    } else { # down
        if(sign(firstslope) < 0 & sign(secondslope) < 0 & abs(ratio) > threshold) {
            return("accelerating")
        } else if(sign(firstslope) > 0 & sign(secondslope) < 0 & abs(ratio) > 0.5*threshold) {
            return("reversing")
        } else {
            return("unsuccessful")
        }
    }
}


### modifications of get_tau() and classify() to work with x = t instead x = cparam

get_idx_time <- function(X, mintime = 5) { # must set a mintime to discard transients
    idx <- switch(
        direction,
        down = which(apply(X, 1, min) > basin & X[, "time"] >= mintime),
        up = which(apply(X, 1, max) < basin &  X[, "time"] >= mintime)
    )

    deltaT <- attr(X, "deltaT")
    idx <- idx[seq(1, length(idx), by = 1/deltaT)]#(1/deltaT)/10)]

    return(idx)
}

get_samples_time <- function(X, which = c("near", "far"), n = 5, mintime = 5) {
    whichslope <- match.arg(which)

    idx <- get_idx_time(X, mintime = mintime)

    switch(
        whichslope,
        near = rev(idx[seq(length(idx), by = -1, length.out = n)]),
        far = idx[seq(1, by = 1, length.out = n)]
    )
}

get_tau_time <- function(X, EWS, mintime = 5) { # no adjusting sign because t always grows and EWS should also grow
    idx <- get_idx_time(X, mintime = mintime)

    times <- X[idx, "time"]
    ews <- EWS[idx]

    tau <- cor(times, ews, method = "kendall")
    return(tau)
}

get_slope_time <- function(X, EWS, which = c("near", "far"), n = 5, mintime = 5, return.model = FALSE) {
    whichslope <- match.arg(which, c("near", "far"))

    samples <- get_samples_time(X, which = whichslope, n = n, mintime = mintime)

    times <- X[samples, "time"]
    ews <- EWS[samples]
    m <- lm(ews ~ times)

    if(return.model) {
        return(m)
    } else {
        return(as.numeric(coef(m)[2]))
    }
}

classify_time <- function(X, EWS, threshold = 2, n = 5, mintime = 5) {
    firstslope <- get_slope_time(X, EWS, "far", n = n, mintime = mintime)
    secondslope <- get_slope_time(X, EWS, "near", n = n, mintime = mintime)

    firstsign <- sign(firstslope)
    secondsign <- sign(secondslope)
    ratio <- secondslope/firstslope

    ## no sense of direction here
    ## direction <- attr(X, "direction")

    ## if(direction == "up") {
    if(sign(firstslope) > 0 & sign(secondslope) > 0 & abs(ratio) > threshold) {
        return("accelerating")
    } else if(sign(firstslope) < 0 & sign(secondslope) > 0 & abs(ratio) > 0.5*threshold) {
        return("reversing")
    } else {
        return("unsuccessful")
    }
    ## } else { # down
    ##     if(sign(firstslope) < 0 & sign(secondslope) < 0 & abs(ratio) > threshold) {
    ##         return("accelerating")
    ##     } else if(sign(firstslope) > 0 & sign(secondslope) < 0 & abs(ratio) > 0.5*threshold) {
    ##         return("reversing")
    ##     } else {
    ##         return("unsuccessful")
    ##     }
    ## }
}

## helper function to properly capitalize Unsuccessful etc., working from toupper
tosentence <- function(x) {
    if(length(x) > 1) stop("Not currently vectorized.")
    if(nchar(x) == 1) {
        toupper(x)
    } else {
        paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
    }
}
