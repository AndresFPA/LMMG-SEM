# LMMG-SEM
This repository contains the files related to the simulation study conducted in the paper "How do Structural Relations Differ Across Groups and Change Over Time? Introducing Longitudinal Mixture Multigroup Structural Equation Modeling for Longitudinal Multigroup Data". The article aimed to integrate MMG-SEM (Perez Alonso et al., 2024) with Latent Markov Models for longitudinal data analysis.

### Simulation
This folder contains the simulation files. The main files are DatGen.R and do_sim_HPC.R, which generate the data and run the simulation, respectively. Note that do_sim_HPC.R is written specifically for the KU Leuven HPC. The files within the "Analyses" folder are the scripts used to analyze the data. 

### hmm-mmgsem
This folder contains the required scripts to run LMMG-SEM. The main function is in the hmm_mmgsem.R file, but it requires all the other secondary functions to work correctly.

### Empirical Example
Includes the files for the re-analysis of the ESS data used by Borgonovi & Pokropek (2019). Note that the data files are not directly included in this repository, but they can be openly accessed at the ESS Data Portal (https://ess.sikt.no/en/).

# References
Borgonovi, F., & Pokropek, A. (2019). Education and Attitudes Toward Migration in a Cross Country Perspective. Frontiers in Psychology, 10. https://doi.org/10.3389/fpsyg.2019.02224

Perez Alonso, A. F., Rosseel, Y., Vermunt, J. K., & De Roover, K. (2024). Mixture multigroup structural equation modeling: A novel method for comparing structural relations across many groups. Psychological Methods. https://doi.org/10.1037/met0000667
