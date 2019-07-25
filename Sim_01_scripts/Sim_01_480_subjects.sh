#!/bin/bash
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l h_rt=23:59:00
#$ -t 1:100
#$ -cwd
#$ -o $HOME/log
#$ -e $HOME/log

. /etc/profile

module add matlab

matlab -nodisplay -nojvm -r Sim_01_480_subjects

