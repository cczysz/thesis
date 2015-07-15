#!/bin/bash

#PBS -N czysz_apt_compare
#PBS -S /bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb

#PBS -o $HOME/mask.out
#PBS -e $HOME/job.err

module load R/3.1.0

R CMD BATCH --no-save --no-restore /home/t.cczysz/thesis/scripts/compare_apt.R
