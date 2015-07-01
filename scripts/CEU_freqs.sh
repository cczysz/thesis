#!/bin/bash

#PBS -N czysz_ceu_freqs
#PBS -S /bin/bash
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=12gb

#PBS -M czysz@uchicago.edu

#PBS -o $HOME/
#PBS -e $HOME/job.err

module load vcftools
module load tabix

# DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
# SCRIPT_DIR=/home/t.cczysz/thesis/scripts
EXON_PROBES=/group/stranger-lab/moliva/ImmVar/probes_mapping/intersections/probes.gencode.v22.exon.intersect.bed

tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz -B $EXON_PROBES > /home/t.cczysz/genotypes.vcf
