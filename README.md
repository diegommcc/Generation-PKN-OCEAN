# Generation of prior knowledge network for OCEAN

This repository contains a set of functions to generate organism-specific prior knowledge networks (PKNs) for the OCEAN R package. It is structured as follows: 

* `data`: IT contains the starting data needed to generate PKNs for every organism. These data come from [this repository](https://github.com/SysBioChalmers/), a platform of genome-scale metabolic models (GEMs) for different organisms. For more details about how these GEMs were generated, see Wang et al., 2021. In this repository, GEMs for human (_Homo sapiens_), mouse (_Mus musculus_) and rat (_Ratus norvegicus_) are available, but the rest can be easily obtained. In addition, information to determine cofactors is also already downloaded and stored as an RDS file.
* `notebooks`: Rmd with a descriptive analysis of the mouse PKN. 
* `output`: final PKNs in form of RDS files. 
* `reports`: HTML and plots (png) generated after render Rmds in `notebooks`. 
* `src`: there are two R files: 
    * `final_functions_PKN_OCEAN.R`: contains the set of functions in charge of generating the PKNs. 
    * `generation_networks_OCEAN.R`: script which loads `final_functions_PKN_OCEAN.R` functions and write the final PKNs as RDS files in `output`. 

These functions will be available through the `OmnipathR` R package in the future. 

# References

* Wang H, Robinson JL, Kocabas P, Gustafsson J, Anton M, Cholley PE, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. Proceedings of the National Academy of Sciences. 2021 Jul 27;118(30):e2102344118.

