#!/bin/bash
#PBS -l nodes=4:ppn=8
#PBS -l walltime=12:00:00
#PBS -q pace-ice
#PBS -N {NAME}
#PBS -o stdout
#PBS -e stderr
#PBS -j oe
#PBS -m abe
#PBS -M luohaozh@gmail.com

cd {SCRIPT_LOCATION}
source /gpfs/pace1/project/chbe-medford/medford-share/envs/espresso-5.1.r11289-$
{COMMAND}