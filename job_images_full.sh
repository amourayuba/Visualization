#!/bin/bash
#SBATCH -t 0-6:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=120G
#SBATCH --nodes=1

source /home/ayuba/projects/def-taylor/ayuba/halo_profiles/bin/activate
python make_images_full_nogroup.py "m25s7" 0 2048 "3.2"



