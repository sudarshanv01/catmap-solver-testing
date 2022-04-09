#!/bin/bash

# CASES=("co_oxidation" "co_oxidation_ads_ads" "co_hydrogenation" "methanol_synthesis")
CASES=("co_oxidation_ads_ads")
SOLVERS=("coverages" "numbers_fix_xstar" "numbers_free_xstar")
BASE=$1
TEMPLATE="input_files"
INPUT_PREFIX="input"
ENERGY_PREFIX="energies"
SOLVER_SP_PREFIX="solver_specifics"
MID="_"

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