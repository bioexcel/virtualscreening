#!/bin/bash
#SBATCH --job-name="pymdsetup"
#SBATCH -D .
#SBATCH --output=serial_%j.out
#SBATCH --error=serial_%j.err
#SBATCH --nodes=3
#SBATCH --time=00:20:00
#SBATCH --qos=bsc_ls
python workflows/gromacs_full.py workflows/conf_2mut_nt0.yaml mare_nostrum 3
