#!/bin/bash

#PBS -N czysz_mod_r
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb

#PBS -o $HOME/
#PBS -e $HOME/job.err

module load R/3.1.0

DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
R CMD BATCH --no-save --no-restore $DATA_DIR/mod_de.r
