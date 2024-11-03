# spatialEWS

This repository contains code and data in support of our manuscript, "Applicability of spatial early warning signals to complex network dynamics." The arXiv version of our manuscript is available [here](https://arxiv.org/abs/2410.04303). Please cite this article when you use this code.

This README is divided into four parts. First, we demonstrate conducting ascending and descending sequences of simulations on a relatively small network. These simulations should be possible on a standard laptop; we tested these simulations on an x86 64-bit laptop with four Intel i3-5010U CPUs at 2.10GHz. In particular, we address finding reasonable parameters for a given network. Second, we show how we compute each of the four "spatial" early warning signals (EWSs) on the data we produce in our simulations. Third, we show how to evaluate each EWS using Kendall's $\tau$ and our classification algorithm. Finally, we describe how to produce rough versions of the figures we included in the manuscript.

## Simulation sequences

To demonstrate our simulation procedure, we use a Barab\'asi-Albert network with 100 nodes. We included this network in this repository, so you can use it. We also need the adjacency matrix of the network and the number of nodes. Note that the steps below paraphrase from `simulate-model.R`, which is intended to conduct simulations as controlled by a shell script to support large simulations on a computing cluster. 

```R
library(igraph)

modelnetworks <- readRDS("./data/modelnetworks.rds")
g <- modelnetworks$barabasialbert
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)
```

We have produced a small R package called [`sdn`](https://github.com/ngmaclaren/sdn) which implements stochastic differential equations on networks and the dynamical models we use in this manuscript. Please follow the instructions in that repository to follow the remaining simulation steps. 

We will rely on parallel processing and our sdn package for the remaining simulation steps.

```R
library(parallel) # NB: we are testing on a Linux machine. If you are using Windows, you may need a different parallelization method.
ncores <- detectCores() - 1
library(sdn)
```

For demonstration purposes, here we use the coupled double-well dynamics. (This dynamics is implemented in sdn.) These dynamics are implemented as functions in the [deSolve](https://cran.r-project.org/package=deSolve) style for interoperability. Please see the sdn documentation.

```R
model <- doublewell
modelparams <- .doublewell
```

Initially, we will use the global stress parameter $u$ as the control parameter and conduct ascending simulations. We will also start with some initial guesses at parameter values. 

```R
cparam <- "u"
lmax <- 100 # how many samples to take from the control parameter range
rng <- seq(0, 5, length.out = lmax) # the range of control parameter values considerd is [0, 5]
D <- 0.05
xinit <- rep(modelparams$xinit.low, N)

## for conducting simulations with sdn:
deltaT <- 0.01
simtime <- 50 # 50 user time units
params <- c(modelparams, list(A = A))
control <- list(ncores = ncores, times = 0:50, deltaT = deltaT)
```

Although the coupled double-well dynamics does not need it, we want to highlight that dynamics which involve an equilibrium at zero (e.g., the SIS dynamics) need special consideration. We tell sdn to set any value below zero to zero through an argument in the control list:

```R
if(args$model %in% c("SIS", "mutualistic", "genereg")) {
    control$absorbing.state <- list(value = 0, which = "floor")
}
```

Now, we can run a test simulation like this

```R
system.time(
    X <- sde(xinit, control$times, model, params, control)
)
```

which took 0.651 seconds on our machine. This `X` is in the output style of deSolve: the first column is the simulation time step and the remaining columns are the $x_i$ values for each node at that time step. Note that this `X` is for one particular control parameter value, namely, $D=0.05$ and $u=0$:
```R
params$u
params$D
```

To simulate across our parameter range, we use sdn's `solve_in_range()`:

```R
system.time(
    result <- solve_in_range(rng, cparam, model, xinit, params, control, kind = "sde")
)
```

That ran in 29.7 seconds on our machine. 

Let's look at what we produced using a diagnostic plot from sdn:

```R
bifplot(result, rng, col = adjustcolor(1, 0.5), lwd = 0.5)
```

While it's nice to see the complete transition, for evaluating EWSs in our manuscript, we are only concerned with the control parameter range up to the first transition of any node. We can update our simulation to focus on this range, thus providing finer resolution near the bifurcation.

```R
rng <- seq(0, 2, length.out = lmax)
result <- solve_in_range(rng, cparam, model, xinit, params, control, kind = "sde")
bifplot(result, rng, col = adjustcolor(1, 0.5), lwd = 0.5)
```

Perhaps $u \in [0, 1.5]$ would be better, but $[0, 2]$ matches what we used for the manuscript (see "Simulation parameters.ods", which is referred to from the Supplementary Materials of the manuscript). Note that this kind of procedure---proposing initial parameter values, checking with some diagnostic, then rerunning the simulations---will need to be done for each new network and for each simulation condition. For large networks, this procedure may be best done by parallelizing across many cores, probably with a high performance computing cluster. 

To see the descending simulations, we update our parameters and re-run:
```R
xinit <- rep(modelparams$xinit.high, N)
rng <- seq(0, -4, length.out = lmax)
result <- solve_in_range(rng, cparam, model, xinit, params, control, kind = "sde")
bifplot(result, rng, col = adjustcolor(1, 0.5), lwd = 0.5)
```

As noted above, we provide `simulate-model.R` to perform these functions from the command line. An example script which produces simulations for each of our ten simulation conditions for one network (i.e., the so-called GKK network) is found in `add-gkk.sh`. For working on the cluster, please see `request-simulation.sh`, which schedules one simulation sequence to be run as one job, and `individual-request.sh` and `build-simulation-request.sh`, which send one and many, respectively, simulation requests to the scheduler. Simulation output from `simulate-model.R` also includes several attributes that assist in follow-on analyses, so using it in an interactive session can be useful.


## Computing early warning signals

Moran's $I$ is defined as 

$$I_{\rm M} = \frac{N}{W} \frac{\sum_{i=1}^N \sum_{j=1}^N A_{ij} (x_i - \overline{x}) (x_j - \overline{x})}{\sum_{i=1}^N (x_i - \overline{x})^2}$$.

Our implementation of Moran's $I$ is called `global_moran()` and is in `./calc-functions.R`.

We use the base R implementation of the sample standard deviation,

$$s = \sqrt{\frac{1}{N-1} \sum_{i=1}^N (x_i - \overline{x})^2}$$.

There are several ways to compute skewness and kurtosis---we use the method based on central moments. The $k$th central moment of a vector of values $\mathbold{x}$ is

$$m_k = \frac{1}{N} \sum_{i=1}^N (x_i - \overline{x})^k$$.

Skewness is based on the third central moment and is defined as

$$g_1 = \frac{m_3}{m_2^{3/2}}$$.

Because we investigate both ascending and descending simulations, for which the expected change in shape of the distribution is opposite (i.e., becoming more skewed as a control parameter increases and become less skewed as a control parameter decreases) we define $$g_1' \equiv g_1$$ for ascending simulations and $$g_1' \equiv -g_1$$ for descending simulations. 

Kurtosis is defined as

$$g_2 = \frac{m_4}{m_2^2}$$. 

We use the implementations of skewness and kurtosis in the [moments](https://cran.r-project.org/package=moments) package. The source code is available at that link: see `./R/skewness.R` and `./R/kurtosis.R` in the directory of that package. 

To compute EWSs, produce a data matrix in which each column corresponds to a node of the network and each row corresponds to a control parameter value. Then, if such a matrix is called `result` (for example, the `result` we computed above), then we run

```R
source("calc-functions.R") # for global_moran()
apply(result, 1, global_moran, A = A)
apply(result, 1, sd)
apply(result, 1, moments::skewness) # not sign-corrected; see below for an example of how to correct the sign
apply(result, 1, moments::kurtosis)
```

This procedure computes one EWS value for each row of the matrix, corresponding to one control parameter value. In our manuscript, we were only interested in EWS values before the first transition of any node. We defined a set of variables and functions to easily identify those values (see `./calc-functions.R`). 

First, we define an $x_i$ value which will serve as the limit of an initial basin of attraction. For example, for our coupled double-well dynamics, we use $(r_1,r_2,r_3)=(1,3,5)$. In the absence of noise and coupling between nodes, $r_2$ is an unstable equilibrium. Specifically, if $r_1 \leq x_i < r_2$, then $x_i$ will move towards $r_1$. If $r_2 < x_i \leq r_3$, then $x_i$ will move towards $r_3$. In the presence of noise and coupling, equilibrium values of $x_i^*$ are not exactly equal to $r_1$ or $r_3$. However, we consider $x_i < r_2$ to be in the lower state and $x_i > r_2$ to be in the upper state. We define a global variable called `basins` in `./calc-functions.R` which stores the equivalent boundary values for each of the four dynamics models we used.

Second, we define a function `get_idx()` which returns the row indices of a specially defined data frame which correspond to the values of the control parameter for which all $x_i$ are in their original state (i.e., are on the same side of the boundary value as their initial value). Our function `get_idx()` relies on several pieces of information which `simulate-model.R` stores alongside the simulation output, including the direction of the simulation sequence, which is relevant for deciding whether or not an $x_i$ value is still in its original basin of attraction. Because we allow for a list of several repetitions of a given simulation sequence, we also include a function `promote_df()` to handle the output of `simulate-model.R`. To compute EWSs for a single simulation sequence, follow these steps:

```R
## Assuming that `simresult` is the output of `simulate-model.R`...
df <- promote_df(simresult)
A <- attr(df, "params")$A # retrieve the adjacency matrix from the "promoted" df

## if using `result` from above in this README, need to add some attributes to satisfy calc-functions.R functions, which assume the output of simulate-model.R:
attr(result, "direction") <- "down" # our most recent `result` started from xinit.high and made u more negative
attr(result, "model") <- "doublewell"
attr(result, "bparam.vals") <- rng
df <- result

## Continuing, same for both methods...
moranI <- apply(df, 1, global_moran, A = A)
ssd <- apply(df, 1, sd)
skew <- apply(df, 1, moments::skewness)
if(attr(df, "direction") == "down") skew <- -skew
kurt <- apply(df, 1, moments::kurtosis)
```

## Evaluate early warning signals

Given the object `df` above and the computed EWSs, computing Kendall's $\tau$ is relatively easy. We provide a function `get_tau()` in `calc-functions.R` which computes the $\tau$ values only for the latter half of the "home range" of the simulations. (See the manuscript for the definition of home range.)

```R
taus <- list(
    moranI = get_tau(df, moranI, adjust.sign = TRUE),
    ssd = get_tau(df, ssd, adjust.sign = TRUE),
    skew = get_tau(df, skew, adjust.sign = TRUE),
    kurt = get_tau(df, kurt, adjust.sign = TRUE)
)
```

For our classification algorithm, we select a small number of observations far from the first bifurcation and a small number of observations near it. We provide functions to do these tasks in `calc-functions.R`: 
- `get_samples()` returns the indices for the first $n$ control parameter values or the last $n$ of them before the transition, depending on arguments.
- `get_slope()` uses `get_samples()` and returns either just the slope of the EWS (with respect to the control parameter) at the far and near points, or the whole linear model as an object.
- `classify()` uses the slopes returned by `get_slope()` and applies our decision rule, returning a character vector (length 1) with the classification result (accelerating, reversing, or unsuccessful).

For example, 
```R
n <- 5
classification <- list(
    moranI = classify(df, moranI, n = n),
    ssd = classify(df, ssd, n = n),
    skew = classify(df, skew, n = n),
    kurt = classify(df, kurt, n = n)
)
```

## Manuscript figures

For copyright purposes, we are not making all networks available in this repository. The networks can be downloaded from their original sources or by contacting me. EWS data for making figures etc. is in `./data/EWS-data.RData`. This `.RData` file contains the EWSs computed from each simulation. The files to produce a raw version of manuscript figures 1, 2, 4, 5, and S1 are provided:
- Figure 1: `example-method.R`
- Figure 2: `tauplot.R`
- Figure 4: `example-drug.R`
- Figure 5: `example-drug.R`
- Figure S1: `example-lattice.R`
Create a subdirectory in your local clone called `./img/` before running those files.

To reproduce `./data/EWS-data.RData` itself (i.e., to recompute the EWS), extract `./data/sims.tar` to a directory called `./data/sims/` (there will be 360 simulation output `.rds` files) and run `./prepare-EWS-data.R`. 

Finally, to re-run our classification results and produce the raw data for Figure 3 and the numerical results presented in the manuscript, section II B, see `./classification-results.R`.

Both `./prepare-EWS-data.R` and `./classification-results.R` take several seconds to run on our machine. Additionally, `./prepare-EWS-data.R` saves the current R environment, whatever it is, as a side effect. We recommend that you run `./prepare-EWS-data.R` in a clean session.
