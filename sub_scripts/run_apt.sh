#!/bin/bash

#PBS -N czysz_apt
#PBS -S /bin/bash
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb

#PBS -o $HOME/apt.out
#PBS -e $HOME/apt_job.err

APT_DIR=/home/t.cczysz/apt/apt-1.17.0-x86_64-intel-linux/bin

$APT_DIR/apt-probeset-summarize -a rma-mask \

