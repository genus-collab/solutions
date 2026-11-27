#!/bin/bash
#SBATCH--comment=PicardFuchs
#SBATCH--nodes=1
#SBATCH--ntasks-per-node=1
#SBATCH--ntasks=1
#SBATCH--cpus-per-task=1
#SBATCH--mem=64G
#SBATCH--output=PicardFuchs.out
#SBATCH--time=96:00:00
module load Mathematica/13.2.1
wolframscript -file genRmat_g4_partA_pf.m 9
