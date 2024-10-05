#!/bin/bash

model="mutualistic"
network="football"
bparam="u"
direction="down"
uinit="-5"
ntrials=1

jobname=${network}-${model}-${bparam}-${direction}

sbatch --job-name=${jobname} --output=./outfiles/${jobname}.out request-simulation.sh ${model} ${network} ${bparam} ${direction} ${uinit} ${ntrials}
