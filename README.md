# spatialEWS

Assessment of spatial early warning signals. 

Note that for copyright purposes we are not making all networks available in this repository---the networks can be downloaded from their original sources or by contacting me. EWS data for making figures etc. is in `./data/EWS-data.RData`. This `.RData` file contains the computed early warning signals from each simulation. The files to produce raw version of manuscript figures 1, 2, 4, 5, and S1 are provided:
- Figure 1: `example-method.R`
- Figure 2: `tauplot.R`
- Figure 4: `example-drug.R`
- Figure 5: `example-drug.R`
- Figure S1: `example-lattice.R`
Create a subdirectory in your local clone called `./img/` before running those files. 

## Computing early warning signals

Our implementation of computing Moran's I is in `calc-functions.R`. To compute early warning signals, produce a data matrix in which each column corresponds to a network node and each row corresponds to a control parameter value. Then, if such a matrix is called `X`,

```R
source("calc-functions.R") # for global_moran()
apply(X, 1, global_moran, A = A)
apply(X, 1, sd)
apply(X, 1, moments::skewness)
apply(X, 1, moments::kurtosis)
```

## Conducting simulations

The base file that performs a sequence of simulations with a given set of parameters is `simulate-model.R`. This file relies on a separate package which can be found [here](https://github.com/ngmaclaren/sdn). 

- `add-gkk.R` is a script that does all ten simulation sequences for one network (in this case, the GKK network) on a local machine using parameters stored in `data/simulation-parameters.csv`. 
- The scripts `individual-request.sh`, `request-simulation.sh`, and `build-simulation-request.sh` are intended for use on a high-performance computing cluster. As written, they are designed for the Slurm workload manager as implemented for a particular computing cluster and may need modification to work in a different environment. 

## Workflow

1. Create a network, e.g.:
```R
library(igraph)
set.seed(123)

N <- 100
M <- 300

g <- largest_component(sample_fitness_pl(N, M, 2))

```
2. Identify relevant parameter ranges.

Run `simulate-model.R` in an interactive session. 
Starting on line 62, update to the data you want, in this case (other options in comments):

```R
if(interactive()) {
    args$network <- "gkk"
    args$model <- "doublewell" # "genereg" "SIS" "doublewell" "mutualistic"
    args$bparam <- "D" # "u"
    args$direction <- "up" # "down"
}

```
## Depends

deSolve, sdn
