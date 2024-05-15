#!/bin/bash

#SBATCH --time=01:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=100000

#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=NONE

module load gcc openmpi r

model=$1
network=$2
bparam=$3
# simkind and simtime can be default
direction=$4
# sigma, xinit, and Dinit can be default
uinit=$5
ntrials=$6

Rscript simulate-model.R  --model=${model} --network=${network} --bparam=${bparam} --direction=${direction} --uinit=${uinit} --ntrials=${ntrials}
