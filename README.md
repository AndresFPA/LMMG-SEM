# LMMG-SEM
This repository contains the files related to the simulation study conducted in the paper "TITLE". The article aimed to integrate MMG-SEM with Latent Markov Models for longitudinal data analysis.

## Simulation
This folder contains the simulation files. The main files are DatGen.R and do_sim_HPC.R, which perform the data generation and simulation, respectively. Note that do_sim_HPC.R is written specifically to work with the HPC from KU Leuven. The files within the "Analyses" folder are the scripts used to analyze the data. 

## hmm-mmgsem
This folder contains the required script to run LMMG-SEM. The main function is in the hmm_mmgsem.R file, but it requires all the other secondary function to work correctly.
