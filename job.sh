#!/bin/bash
#PBS -l ncpus=1
#PBS -l mem=64GB
#PBS -l jobfs=180GB
#PBS -q normal
#PBS -P zv32
#PBS -l walltime=36:00:00
#PBS -l storage=gdata/zv32+scratch/zv32
#PBS -M Quanling.Deng@anu.edu.au
#PBS -l wd

# Now list your executable command (or a string of them).                                                                
# Example for code compiled with a software module:                                                                      

module load matlab/R2020b
module load matlab_licence/anu
matlab -nodisplay -nosplash -r "outputDir='$PBS_JOBFS',numberOfWorkers=$PBS_NCPUS, EuDAloop3nx, quit" > $\
PBS_JOBID.log

