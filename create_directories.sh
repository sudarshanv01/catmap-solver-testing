#!/bin/bash

CASES=("co_oxidation" "co_oxidation_ads_ads" "co_hydrogenation" "methanol_synthesis" "co_reduction")
SOLVERS=("coverages" "numbers_free_xstar")
TEMPLATE="input_files"
INPUT_PREFIX="input"
ENERGY_PREFIX="energies"
SOLVER_SP_PREFIX="solver_specifics"
MID="_"

if [ $# -eq 0 ]
  then
    echo "No arguments supplied; please supply the convergence of the solver name."
    exit 0
fi

BASE=$1

if [ "$BASE" == "tex" ]
	then
		echo "Creating directories for text output."
		for CASE in ${CASES[@]}; do
			for SOLVER in ${SOLVERS[@]}; do
				mkdir -p $BASE/$CASE/$SOLVER
				# Copy input files
				cp $TEMPLATE/$INPUT_PREFIX$MID$CASE.mkm $BASE/$CASE/$SOLVER/$INPUT_PREFIX.mkm 
				cp $TEMPLATE/$ENERGY_PREFIX$MID$CASE.txt $BASE/$CASE/$SOLVER/$ENERGY_PREFIX.txt 
				# Copy python job files
				cp $TEMPLATE/tex_mkm.py $BASE/$CASE/$SOLVER/tex_mkm.py
				# Copy the json file required for the solver
				cp $TEMPLATE/$SOLVER_SP_PREFIX$MID"low".json $BASE/$CASE/$SOLVER/$SOLVER_SP_PREFIX.json
			done
		done
else
	for CASE in ${CASES[@]}; do
		for SOLVER in ${SOLVERS[@]}; do
		mkdir -p $BASE/$CASE/$SOLVER
		# Copy input files
		cp $TEMPLATE/$INPUT_PREFIX$MID$CASE.mkm $BASE/$CASE/$SOLVER/$INPUT_PREFIX.mkm 
		cp $TEMPLATE/$ENERGY_PREFIX$MID$CASE.txt $BASE/$CASE/$SOLVER/$ENERGY_PREFIX.txt 
		# Copy python job files
		cp $TEMPLATE/mkm_job.py $BASE/$CASE/$SOLVER/mkm_job.py
		# Copy the json file required for the solver
		cp $TEMPLATE/$SOLVER_SP_PREFIX$MID$BASE.json $BASE/$CASE/$SOLVER/$SOLVER_SP_PREFIX.json
		done
	done
fi
