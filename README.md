# Visualization Simulations of large scale structures 

# Description 
This is a set of codes aimed at visualizing N-body simulations. These are a large number of points, arranged 
in multi-layered structures, and require specific smoothing codes to illustrate the features and compute the 
density fields. 

There is also a set of codes to make a video of the evolution of the large scale structures with time. And for 
three levels of zooms. 

# Requirements 
Python.
Because the following scripts are run through cloud computing, all the slrum scripts 
used in this project assume the existence of a python virtual environnement that is loaded 
through: 
source $ENV/bin/activate 
where $ENV is where the environnement is located. This python environnement has to have the following 
modules isntalled and active:  
readgadget and MAS_library from Pylians https://pylians3.readthedocs.io/en/master/
scipy, numpy, matplotlib, pandas, astropy
From vieo making: 
opencv 

# Data requirements 
Gadget snapshots: in basefolder/sim/snapdir_**
Halo files (if needed) in basefolder/sim/AHF/halos/
Basefolder variable can be changed python script. 
The scripts illustrating halo positions on top of density fields
(make_images_full.py) will also require [sim]_prefixes.txt file with the name 
of the AHF prefixes of each snapshot in basefolder/sim/. 

# How to use 
The functions needed to read, smooth and visualize are in 
## base_for_images.py
This contains routines to read particle positions, halo positions, calculating densities and 
smoothing fields. 

## How to run the programs to generate the images
Because the number of particles can be very large (1024^3). It will require running 
through cloud computing. There are very short slurm scripts for each zoom level. 

### make_images_full.py/make_images_full_nogroup.py
These codes can be run through slurm script **job_images_full.sh** 
They will require 5 arguments: 
sim: simulation name 
smin: Snapshot where to start making images
grid: There will be grid^3 number of grid cells for smoothing 
vm: type of maximum color density range. If "variable" as the density changes, the color range changes 
with average density. Otherwise, color range is fixed for all snapshots. 
mlim: Halo mass limit to illustrate if put halos.

These will generate images with or without halos above mlim for all snapshot starting from smin.
If one does not want halo positions to be highlighted in the snapshot images, 
need to modify the **job_images_full.sh** sbatch script to put make_images_full_nogroup.py 
instead of make_images_full.py, and remove the fitfh argument mlim. 

### make_images_50Mpc.py 
This code can be run through the slurm batch script job_image_zoom50.sh. It is 
10x zoom-in on a portion of the simulation. It requires four arguments 
sim: simulation name 
smin: Snapshot where to start making images
grid: There will be grid^3 number of grid cells for smoothing 
vm: type of maximum color density range. If "variable" as the density changes, the color range changes 
with average density. Otherwise, color range is fixed for all snapshots. 

### make_images_10Mpc.py 
This code can be run through the slurm batch script job_image_zoom50.sh. It is 
50x zoom-in on a portion of the simulation. It requires two arguments 
sim: simulation name 
smin: Snapshot where to start making images 

There is less optionality for changing smoothing and grids because this level of zoom requires more 
specific fine-tuning optimization. 

## Make videos out the images 
make_videos.ipynb is a short notebook that makes a video out of the images 
produces from different scripts. 

## Tests and checks 
The notebooks testing_visualisation_scripts.ipynb, visualisation_movie_halos_ipynb and 
visualisation_movie_zoom_cluster.ipynb are for making individual images and testing the codes and choosing the colors 
and fine-tuning. 