#!/bin/bash
#$ -S /bin/bash
#$ -l h_vmem=7G
#$ -l h_rt=11:59:00
#$ -t 1:600
#$ -cwd
#$ -o $HOME/log
#$ -e $HOME/log

. /etc/profile

module add matlab

matlab -nodisplay -nojvm -r Sim_05_120_subjects

