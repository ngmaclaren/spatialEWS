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

check_zerocrossing <- function() {}
q1 <- function() {}

check_magnitude <- function() {}
q2 <- function() {}
