#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --partition=short
#SBATCH --mem=90G
#SBATCH --cores-per-socket=20
#SBATCH --output=%A-%a.out
#SBATCH --error=%A-%a.err

module load R/4.4.1

# above specifications used to run all RETRO_scoring_<>.R files in RETRO-analysis repo

Rscript RETRO_scoring_cyc.R
Rscript RETRO_scoring_myo.R
Rscript RETRO_scoring_myelo.R
