#!/bin/bash

#PBS -N czysz_cd14_peer
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb

#PBS -o $HOME/
#PBS -e $HOME/job.err

module load R/3.1.0

DATA_DIR=/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14
R CMD BATCH --no-save --no-restore $DATA_DIR/cd14_peer.r
