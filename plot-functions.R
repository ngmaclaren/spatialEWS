bifplot <- function(X, cparam, showmean = FALSE, ...) {
                                        # X is the output of solve_in_range()
                                        # cparam should be a sequence
    stopifnot(length(cparam) == nrow(X))

    with(list(...), {
        matplot(cparam, X, type = "l", lty = 1, xlab = "Control parameter", ylab = "x", ...)
    })
    if(showmean) lines(cparam, rowMeans(X), col = 2, lwd = 2)
}
