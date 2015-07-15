#!/bin/bash

#PBS -N czysz_apt
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb

#PBS -o $HOME/apt.out
#PBS -e $HOME/apt.err

APT_DIR=/home/t.cczysz/apt/apt-1.17.0-x86_64-intel-linux/bin
AFFY_FILE_DIR=/group/stranger-lab/forCharles/probes_mapping/Hugene_info

HUGENE=HuGene-1_0-st-v1.r4

$APT_DIR/apt-probeset-summarize \
-a rma \
-p $AFFY_FILE_DIR/$HUGENE.pgf \
-c $AFFY_FILE_DIR/$HUGENE.clf \
-b $AFFY_FILE_DIR/$HUGENE.bgp \
--kill-list /home/t.cczysz/probes_to_kill.txt \
-o apt_out \
--cel-files /home/t.cczysz/cd4_celfiles.txt
#  -s /home/t.cczysz/probesets \
# -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish
