#!/bin/bash

#PBS -N czysz_qc_r
#PBS -S /bin/bash
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb

#PBS -o $HOME/
#PBS -e $HOME/job.err

module load R/3.1.0

DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
R CMD BATCH $DATA_DIR/analyze_affy.r
