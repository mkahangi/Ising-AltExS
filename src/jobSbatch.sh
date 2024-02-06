#!/bin/bash
#SBATCH --job-name=Ising_AltExS
#SBATCH --output=OUT_%x
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --gpus=1
#SBATCH --mail-user=mkahangi@central.uh.edu   # Replace with your email address
#SBATCH --mail-type=END                     # Receive an email when the job finishes

make clean
make run
