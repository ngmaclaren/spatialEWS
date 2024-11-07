load("./data/EWS-data.RData")

investigate_lattice <- FALSE

tdat <- data.frame(
    tau.I = taus$moranI, # these are 'sign-adjusted'
    tau.sd = taus$ssd,
    tau.var = taus$var,
    tau.cv = taus$cv,
    tau.sk = taus$skew,
    tau.ku = taus$kurt,
    dynamics = dyns,
    network = nets,
    cparam = bparams,
    direction = directions,
    row.names = seq_along(dfs)
)
tdat$dynamics <- factor(tdat$dynamics, levels = rev(c("doublewell", "mutualistic", "SIS", "genereg")))
tdat$network <- factor(tdat$network)
tdat$cparam <- factor(tdat$cparam)
tdat$direction <- factor(tdat$direction)

tdat <- reshape(
    tdat, varying = c("tau.I", "tau.sd", "tau.var", "tau.cv", "tau.sk", "tau.ku"), # , "tau.sd"
    v.names = "tau", timevar = "EWS", times = c("I", "ssd", "var", "cv", "g1", "g2"), # , "s"
    direction = "long", new.row.names = 1:10000
)

agg <- aggregate(tau ~ dynamics + cparam + direction + EWS, data = tdat, FUN = mean)
reshape(agg, timevar = "EWS", idvar = c("dynamics", "cparam", "direction"), direction = "wide")

agg2 <- aggregate(tau ~ dynamics + cparam + direction + EWS, data = tdat, FUN = function(x) sum(x > 0.7)/length(x))
reshape(agg2, timevar = "EWS", idvar = c("dynamics", "cparam", "direction"), direction = "wide")
