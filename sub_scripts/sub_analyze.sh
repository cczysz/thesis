#!/bin/bash

#PBS -N czysz_new
#PBS -S /bin/bash
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=12gb

#PBS -M czysz@uchicago.edu

#PBS -o $HOME/
#PBS -e $HOME/job.err

module load R/3.1.0

# DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
SCRIPT_DIR=/home/t.cczysz/thesis/scripts
R CMD BATCH $SCRIPT_DIR/cd4_analysis.r
