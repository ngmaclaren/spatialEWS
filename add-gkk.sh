#!/bin/bash

models=("doublewell" "mutualistic" "genereg" "SIS")
network="gkk"
bparams=("D" "u")
directions=("up" "down")

ntrials=1

for model in ${models[@]}; do
    for bparam in ${bparams[@]}; do
	for direction in ${directions[@]}; do
	    if [[ ${model} == "mutualistic" && ${direction} == "up" ]]; then
		continue
	    fi

	    if [[ ${model} == "genereg" && ${direction} == "up" ]]; then
		continue
	    fi

	    if [[ ${model} == "SIS" && ${bparam} == "u" ]]; then
		continue
	    fi

	    if [[ (${model} == "doublewell" || ${model} == "mutualistic") &&
		      ${direction} == "down" && ${bparam} == "D" ]]; then
		uinit="-5"; else
		uinit="0"
	    fi

	    Rscript simulate-model.R --model=${model} --network=${network} --bparam=${bparam} --direction=${direction} --uinit=${uinit} --ntrials=${ntrials}
	done
    done
done

	
