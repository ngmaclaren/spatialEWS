#!/bin/bash

# Need to rerun
# doublewell-ug_village-u-up
# doublewell-wiki_ht-u-up
# genereg-london_transport-D-down
# genereg-SITC-D-down
# genereg-yeast-D-down
# genereg-football-u-down
# genereg-genefusion-u-down
# genereg-jung_c-u-down
# genereg-london_transport-u-down
# genereg-yeast-u-down
# mutualistic-football-u-down

model="mutualistic"
network="football"
bparam="u"
direction="down"
uinit="-5"
ntrials=1

jobname=${network}-${model}-${bparam}-${direction}

sbatch --job-name=${jobname} --output=./outfiles/${jobname}.out request-simulation.sh ${model} ${network} ${bparam} ${direction} ${uinit} ${ntrials}
