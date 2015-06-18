#!/bin/bash

### Set the name of the job, where jobname is a unique name for your job
#PBS -N czysz_meqtl

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=48:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=1:ppn=6

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=40gb

### Set the destination for your program's output.
#PBS -o ~/meqtl.out
#PBS -e ~/meqtl.err

# Load the approprite applications
module load R/3.1.0

# Execute the program
R CMD BATCH --no-save --no-restore /group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4/mateqtl.r
