#!/bin/bash

#PBS -N czysz_all_mask
#PBS -S /bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb

#PBS -o $HOME/mask.out
#PBS -e $HOME/job.err

module load R/3.1.0

#DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
R CMD BATCH --no-save --no-restore /home/t.cczysz/thesis/scripts/cd4_cd14_cauAffy.R
