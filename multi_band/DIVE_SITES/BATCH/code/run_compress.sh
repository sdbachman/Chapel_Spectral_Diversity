#!/bin/bash
#PBS -N compress
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=20G
#PBS -l walltime=02:00:00
#PBS -A NCGD0011
#PBS -q casper@casper-pbs
#PBS -j oe

module load python/2.7.14
ncar_pylib

REPL1
##mpiexec -n 1 ./compress.py -f yay.nc
