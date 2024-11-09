#!/bin/bash
#
#SBATCH --job-name=orthoGran
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --partition prod
#SBATCH --mem 40G
#
source ~/miniconda3/bin/activate geo
python viirs-orthoGranule.py
