global_moran <- function(x, A) {
    N <- length(x)
    stopifnot(N == nrow(A))
    stopifnot(nrow(A) == ncol(A))

    W <- sum(A)
    x <- x - mean(x) # center, for convenience below

    (N/W)*sum(A*outer(x, x))/sum(x^2)
}

                                        # We do not currently use this function.
                                        # Leaving in place in case we decide to do something with local Moran.
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

                                        # This function returns the row indices of the relevant range.
                                        # It implements a decision rule about states which says "I'm in my original
                                        # state unless I am as or more extreme than a certain value."
get_idx <- function(df) {
    direction <- attr(df, "direction")
    dynamics <- attr(df, "model")
    basin <- basins[[dynamics]]

    switch(
        direction,
                                        # These are correct: if the test is "as or more extreme than" then the
                                        # basin limit is the first value in the new regime. This idx is the indices
                                        # of all cparam vals in the orignial regime. 
        down = which(apply(df, 1, min) > basin),
        up = which(apply(df, 1, max) < basin) 
    )
}

                                        # Returns the row indices of the first n (far) or last n (near) samples in
                                        # the relevant range.
get_samples <- function(df, which = c("near", "far"), n = 3) {
    whichslope <- match.arg(which)

    idx <- get_idx(df)

    switch(
        whichslope,
        near = rev(seq(length(idx), by = -1, length.out = n)),
        far = seq(1, by = 1, length.out = n)
    )
}

                                        # returns the Kendall's tau between the control parameter and the EWS.
                                        # The EWS should be computed beforehand. If up, then τ, else τ' (= -τ).
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

                                        # Returns the simple slope for the first/last n data in the relevant range.
                                        # Independent: control parameter. Dependent: EWS. 
                                        # Optionally returns the whole model.
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

                                        # Classifies an EWS as accelerating, reversing, or unsuccessful based on
                                        # the slope of the EWS vs. the control parameter and a threshold value
                                        # (default = 2). 
classify <- function(df, EWS, threshold = 2, n = 5) {
    firstslope <- get_slope(df, EWS, "far", n = n)
    secondslope <- get_slope(df, EWS, "near", n = n)

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
