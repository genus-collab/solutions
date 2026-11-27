#!/bin/bash
#SBATCH--comment=GraphSampling
#SBATCH--array=1-50
#SBATCH--nodes=1
#SBATCH--ntasks=1
#SBATCH--cpus-per-task=2
#SBATCH--mem=16G
#SBATCH--output=GraphSampling.out
#SBATCH--time=96:00:00
module load Mathematica/13.2.1
wolframscript -file genRmat_g4_sample.m 4 9 $SLURM_ARRAY_TASK_ID 50 8
