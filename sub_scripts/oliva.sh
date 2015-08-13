#!/bin/bash

#PBS -N czysz_mod_oliva
#PBS -S /bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=16gb

#PBS -o $HOME/
#PBS -e $HOME/mod_oliva.err

module load R/3.1.0

R CMD BATCH --no-save --no-restore /home/t.cri.cczysz/thesis/scripts/generateFilteredGeneExp.R
