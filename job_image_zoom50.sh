#!/bin/bash
#SBATCH -t 0-10:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=220G
#SBATCH --nodes=1

source /home/ayuba/projects/def-taylor/ayuba/halo_profiles/bin/activate
python make_images_50Mpc.py "m25s7" 40 4048 "7.0"


