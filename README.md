# spatialEWS

This repository contains code and data in support of our manuscript, "Applicability of spatial early warning signals to complex network dynamics." The arXiv version of our manuscript should post next week (week of Oct 6). The repository is a work in progress: please be patient.

This README is divided into four parts. First, we demonstrate conducting ascending and descending sequences of simulations on a relatively small network. These simulations should be possible on a standard laptop; we tested these simulations on an x86 64-bit laptop with four Intel i3-5010U CPUs at 2.10GHz. In particular, we address finding reasonable parameters for a given network. Second, we show how we compute each of the four early warning signals (EWSs) on the data we produce in our simulations. Third, we show how to evaluate each EWS using Kendall's $\tau$ and our classification algorithm. Finally, we describe how to produce rough versions of the figures we included in the manuscript.

## Simulation sequences

To demonstrate our simulation procedure, we use a Barab\'asi-Albert network with 100 nodes. We included this network in this repository, so we can access it easily. We will also need the adjacency matrix and the number of nodes readily available. Note that the steps below paraphrase from `simulate-model.R`, which is intended to conduct simulations as controlled by a shell script to support large simulations on a computing cluster. 

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
library(parallel)
ncores <- detectCores() - 1
library(sdn)
```

For demonstration purposes, we will use the coupled double-well dynamics (this dynamics is implemented in sdn). These dynamics are implemented as functions in the [deSolve](https://cran.r-project.org/package=deSolve) style for interoperability. Please see the sdn documentation.

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

Although the coupled double-well dynamics does not need it, we want to highlight that dynamics which involve a fixed point at zero (e.g., the SIS dynamics) need special consideration. We tell sdn to set any value below zero to zero through an argument in the control list:

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

which took 0.651 seconds on our machine. This `X` is in the output style of deSolve: the first column is the simulation time step and the remaining columns are the $x_i$ values for each node at that time step. Note that this `X` is for one particular control parameter value, namely $D=0.05$ and $u=0$:
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

While it's nice to see the complete transition, for evaluating EWS in this manuscript we are only concerned with the control parameter range up to the first transition of any node. We can update our simulation to focus on this range, thus providing finer resolution near the bifurcation.

```R
rng <- seq(0, 2, length.out = lmax)
result <- solve_in_range(rng, cparam, model, xinit, params, control, kind = "sde")
bifplot(result, rng, col = adjustcolor(1, 0.5), lwd = 0.5)
```

Perhaps $u \in [0, 1.5]$ would be better, but $[0, 2]$ matches what we used for the manuscript (see "Simulation parameters.ods"). Note that this kind of procedure---proposing initial parameter values, checking with some diagnostic, then rerunning the simulations---will need to be done for each new network and for each simulation condition. For large networks, this procedure may be best done parallelizing across many cores, such as is possible on a high performance computing cluster. 

To see the descending simulations, we update our parameters and re-run:
```R
xinit <- rep(modelparams$xinit.high, N)
rng <- seq(0, -4, length.out = lmax)
result <- solve_in_range(rng, cparam, model, xinit, params, control, kind = "sde")
bifplot(result, rng, col = adjustcolor(1, 0.5), lwd = 0.5)
```

As noted above, we provide `simulate-model.R` to perform these functions from the command line. An example script which produces simulations for each of our ten simulation conditions for one network (the GKK network) is found in `add-gkk.sh`. For working on the cluster, please see `request-simulation.sh`, which schedules one simulation sequence to be run as one job, and `individual-request.sh` and `build-simulation-request.sh` which send one and many, respectively, simulation requests to the scheduler. 


## Computing early warning signals

Moran's $I$ is defined as 

$$I_{\rm M} = \frac{N}{W} \frac{\sum_{i=1}^N \sum_{j=1}^N A_{ij} (x_i - \overline{x}) (x_j - \overline{x})}{\sum_{i=1}^N (x_i - \overline{x})^2}$$.

Our implementation of Moran's I is called `global_moran()` and is in `./calc-functions.R`.

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

To compute early warning signals, produce a data matrix in which each column corresponds to a network node and each row corresponds to a control parameter value. Then, if such a matrix is called `result` (for example, the `result` we computed above),

```R
source("calc-functions.R") # for global_moran()
apply(result, 1, global_moran, A = A)
apply(result, 1, sd)
apply(result, 1, moments::skewness)
apply(result, 1, moments::kurtosis)
```

This procedure computes one EWS value for each row of the matrix, corresponding to one control parameter value. For this project we are only interested in EWS values before the first transition of any node. We defined a set of variables and functions to easily identify those values (see `./calc-functions.R`). 

First, we define an $x_i$ value which will serve as the limit of an initial basin of attraction. For example, for our coupled double-well dynamics we use $(r_1,r_2,r_3)=(1,3,5)$. In the absence of noise, $r_2$ is an unstable fixed point: if $r_1 \leq x_i < r_2$ $x_i$ will move towards $r_1$ and $r_2 < x_i \leq r_3$ $x_i$ will move towards $r_3$. In the presence of noise and coupling, equilibrium values of $x_i^*$ will not be exactly equal to $r_1$ or $r_3$. However, we consider $x_i < r_2$ to be in the lower state and $x_i > r_2$ to be in the upper state. We define a global variable called `basins` in `./calc-functions.R` which stores the equivalent boundary values for each of our four dynamics. 

Second, we define a function `get_idx()` which returns the row indices of a specially defined data frame which correspond to the values of the control parameter for which all $x_i$ are in their original state (i.e., are on the same side of the boundary value as their initial value). Our function `get_idx()` relies on several pieces of information which `simulate-model.R` stores alongside the simulation output, including the direction of the simulation sequence, which is relevant for deciding whether or not a an $x_i$ value is still in its original basin of attraction. 

## Evaluate early warning signals

Text about our classification algorithm here.

## Manuscript figures

Note that for copyright purposes we are not making all networks available in this repository. ---the networks can be downloaded from their original sources or by contacting me. EWS data for making figures etc. is in `./data/EWS-data.RData`. This `.RData` file contains the computed early warning signals from each simulation. The files to produce raw version of manuscript figures 1, 2, 4, 5, and S1 are provided:
- Figure 1: `example-method.R`
- Figure 2: `tauplot.R`
- Figure 4: `example-drug.R`
- Figure 5: `example-drug.R`
- Figure S1: `example-lattice.R`
Create a subdirectory in your local clone called `./img/` before running those files.

To reproduce `./data/EWS-data.RData` itself (i.e., to recompute the EWS), extract `./data/sims.tar` to a directory called `./data/sims/` (there will be 360 simulation output `.rds` files) and run `./prepare-EWS-data.R`. 

Finally, to re-run our classification results and produce the raw data for Figure 3 and the numeric results presented in the manuscript, section II B, see `./classification-results.R`.

Note that both `./prepare-EWS-data.R` and `./classification-results.R` take several seconds to run. Additionally, `./prepare-EWS-data.R` saves the current R environment, whatever it is, as a side effect. We recommend that you run `./prepare-EWS-data.R` in a clean session.
