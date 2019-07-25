#!/bin/bash
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l h_rt=04:59:00
#$ -t 1:100
#$ -cwd
#$ -o $HOME/log
#$ -e $HOME/log

. /etc/profile

module add matlab

matlab -nodisplay -nojvm -r Sim_04_60_subjects

