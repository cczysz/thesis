#!/bin/bash

#PBS -N snap_results
#PBS -S /bin/bash
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb

#PBS -M czysz@uchicago.edu

#PBS -o $HOME/
#PBS -e $HOME/job.err

module load python/3.4.3

# DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
SCRIPT_DIR=/home/t.cczysz/thesis/scripts
python $SCRIPT_DIR/getBroadLD.py
