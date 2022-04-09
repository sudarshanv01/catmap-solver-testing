#!/bin/bash
#SBATCH -J job
#SBATCH -p xeon16   ### xeon8; xeon16; xeon24
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 168:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e err

for A in $(seq -w 1 16); 
do
  DIR=$(sed -n "${A}p" run_submit)
  cd $DIR
  python mkm_job.py
  cd - 
done

wait