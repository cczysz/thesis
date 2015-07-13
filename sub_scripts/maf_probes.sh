#!/bin/bash

#PBS -N czysz_maf_overlap
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb

#PBS -o $HOME/maf.out
#PBS -e $HOME/maf.err


WD=/group/stranger-lab/forCharles/probes_mapping/intersections
dbsnp142.overlappingprobes_MAF0.1AFR._EAS_EUR.txt
probes.dbsnp142.intersect.bed

while read line
do
	grep -P "$line\s" $WD/probes.dbsnp142.intersect.bed >> /home/t.cczysz/maf.probes
done <$WD/dbsnp142.overlappingprobes_MAF0.1AFR._EAS_EUR.txt	
