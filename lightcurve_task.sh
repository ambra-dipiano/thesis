#!/bin/bash

#SBATCH -N1
#SBATCH --ntasks-per-node=1
#SBATCH --account=pianoambra@morgana

python /mnt/nvme0n1p1/piano_analysis/working-dir/script_RTAdetection_lightcurve.py

