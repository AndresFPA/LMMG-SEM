# LMMG-SEM
This repository contains the files related to the simulation study conducted in the paper "How do Structural Relations Differ Across Groups and Change Over Time? Introducing Longitudinal Mixture Multigroup Structural Equation Modeling for Longitudinal Multigroup Data". The article aimed to integrate MMG-SEM (Perez Alonso et al., 2024) with Latent Markov Models for longitudinal data analysis.

## Simulation
This folder contains the simulation files. The main files are DatGen.R and do_sim_HPC.R, which generate the data and run the simulation, respectively. Note that do_sim_HPC.R is written specifically for the KU Leuven HPC. The files within the "Analyses" folder are the scripts used to analyze the data. 

## hmm-mmgsem
This folder contains the required scripts to run LMMG-SEM. The main function is in the hmm_mmgsem.R file, but it requires all the other secondary functions to work correctly.

# References
Perez Alonso, A. F., Rosseel, Y., Vermunt, J. K., & De Roover, K. (2024). Mixture multigroup structural equation modeling: A novel method for comparing structural relations across many groups. Psychological Methods. https://doi.org/10.1037/met0000667
