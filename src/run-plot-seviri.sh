#!/bin/bash
#
#SBATCH --job-name=plotSeviri
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --partition prod
#SBATCH --mem 40G
#
source ~/miniconda3/bin/activate geo
python seviri-plotFire.py pdV
python seviri-plotFire.py landes
