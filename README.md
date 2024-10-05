# early-warnings-spatial

Assessment of spatial early warning signals. 

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
