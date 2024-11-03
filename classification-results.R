load("./data/EWS-data.RData")

remove_lattice <- TRUE

results <- data.frame(
    result = c(
        mapply(classify, dfs, moranIs, n = n),
        mapply(classify, dfs, ssds, n = n),
        mapply(classify, dfs, vars, n = n),
        mapply(classify, dfs, cvs, n = n),
        mapply(classify, dfs, skews, n = n),
        mapply(classify, dfs, kurts, n = n)
    ),
    EWS = c(
        rep("moranI", length(dfs)),
        rep("ssd", length(dfs)),
        rep("var", length(dfs)),
        rep("cv", length(dfs)),
        rep("skew", length(dfs)),
        rep("kurt", length(dfs))
    ),
    dynamics = dyns,
    network = nets,
    bparam = bparams,
    direction = directions
)

if(remove_lattice) results <- subset(results, network != "lattice")

print("Moran's I")
with(list(df = subset(results, EWS == "moranI")), table(df$dynamics, df$bparam, df$result, df$direction))
print("Spatial standard deviation")
with(list(df = subset(results, EWS == "ssd")), table(df$dynamics, df$bparam, df$result, df$direction))
print("Variance (moment)")
with(list(df = subset(results, EWS == "var")), table(df$dynamics, df$bparam, df$result, df$direction))
print("Coefficient of variation")
with(list(df = subset(results, EWS == "cv")), table(df$dynamics, df$bparam, df$result, df$direction))
print("Spatial skew")
with(list(df = subset(results, EWS == "skew")), table(df$dynamics, df$bparam, df$result, df$direction))
print("Spatial kurtosis")
with(list(df = subset(results, EWS == "kurt")), table(df$dynamics, df$bparam, df$result, df$direction))


### Overall, skewness is successful x% of the 350 cases, etc. Also give the "winners" of the 10 categories (dynamics/control parameter) (how many times is EWS the winner out of 10). Also rank, say 3, 2, 1, 0 for 1st, 2nd, 3rd, 4th place. Give total rank.

success <- aggregate(
    result ~ EWS + dynamics + bparam + direction,
    data = results,
    FUN = function(x) sum(x != "unsuccessful")/length(x)
)
success$EWS <- factor(success$EWS)
success$dynamics <- factor(success$dynamics)
success$bparam <- factor(success$bparam)
success$direction <- factor(success$direction)

ssplit <- split(success, ~ dynamics + bparam + direction)
ssplit <- lapply(ssplit, function(x) {
    x$rank <- rank(x$result)
    x
})
ranks <- do.call(rbind, ssplit)
rownames(ranks) <- seq(nrow(ranks))

ranks$rank <- ranks$rank - 1
                                        # wins
aggregate(rank ~ EWS, data = ranks, function(x) sum(x == 5))
tapply(ranks, ~ EWS + dynamics + bparam + direction, FUN = function(x) max(x$rank, na.rm = TRUE))
                                        # total rank
aggregate(rank ~ EWS, data = ranks, sum)
                                        # overall percent
aggregate(result ~ EWS, data = results, FUN = function(x) sum(x != "unsuccessful")/length(x))
