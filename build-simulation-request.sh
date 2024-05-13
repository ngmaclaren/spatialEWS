#!/bin/bash

# debugging
# models=("doublewell")
# networks=("dolphin")
# bparams=("D")
# directions=("up")

# lists of conditions
# models=("doublewell" "mutualistic" "genereg" "SIS")
# models=("doublewell" "mutualistic" "SIS")
models=("genereg")
networks=(`cat ./data/networknames.txt`)
bparams=("D" "u")
directions=("up" "down")

# standard parameter(s)
ntrials=1 # 50

# for distinguishing job names
counter=1

# loop through all combinations of conditions
for network in ${networks[@]}; do
    # doublewell: all combinations
    for model in ${models[@]}; do
	for bparam in ${bparams[@]}; do
	    for direction in ${directions[@]}; do
		# mutualistic: down/u, down/D
		if [[ ${model} == "mutualistic" && ${direction} == "up" ]]; then
		    continue
		fi

		# genereg: down/u, down/D
		if [[ ${model} == "genereg" && ${direction} == "up" ]]; then
		    continue
		fi

		# SIS: up/D, down/D
		if [[ ${model} == "SIS" && ${bparam} == "u" ]]; then
		    continue
		fi

		# negative u required for doublewell/mutualistic for down/D
		if [[ (${model} == "doublewell" || ${model} == "mutualistic") &&
			  ${direction} == "down" && ${bparam} == "D" ]]; then
		    uinit="-5"; else
		    uinit="0"
		fi

		# jobname=$(printf "sim-%03d" "${counter}")
		jobname=${network}-${model}-${bparam}-${direction}
		# echo ${jobname}
		# echo ${counter}

		sbatch --job-name=${jobname} --output=./outfiles/${jobname}.out request-simulation.sh ${model} ${network} ${bparam} ${direction} ${uinit} ${ntrials}

		((counter++))
	    done
	done
    done
done




# for model in ${models[@]}; do
#     for network in ${networks[@]}; do
# 	for bparam in ${bparams[@]}; do
# 	    for direction in ${directions[@]}; do
# 		# Skip conditions that aren't physical
# 		# if [[ ${model} == "SIS" && ${direction} == "down" ]]; then
# 		#     continue
# 		# fi
# 		if [[ ${model} == "SIS" && ${bparam} == "u" ]]; then
# 		    continue
# 		fi
# 		# if [[ ${model} == "genereg" && ${direction} == "up" ]]; then
# 		#     continue
# 		# fi

# 		# Accommodate negative u required for doublewell/mutualistic models to bifurcate when the bparam is D and the direction is down
# 		if [[ (${model} == "doublewell" || ${model} == "mutualistic") &&
# 			  ${direction} == "down" && ${bparam} == "D" ]]; then
# 		    uinit="-5"; else
# 		    uinit="0"
# 		fi
		    
# 		jobname=$(printf "sim-%03d" "${counter}")

# 		sbatch --job-name=${jobname} --output=./outfiles/${jobname}.out request-simulation.sh ${model} ${network} ${bparam} ${direction} ${uinit} ${ntrials}

# 		((counter++))
# 	    done
# 	done
#     done
# done

