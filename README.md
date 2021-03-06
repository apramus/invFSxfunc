# Data and Code from Ramus et al. (2017) *PNAS*

This repository contains the plot-level data used to generate all figures and code to replicate the primary analysis of ecosystem multifunctionality as presented in

Ramus AP, Silliman BR, Thomsen MS, Long ZT (2017) An invasive foundation species enhances multifunctionality in a coastal ecosystem. *Proc Natl Acad Sci USA* 114(32):8580-8585. http://www.pnas.org/content/114/32/8580

`The authors request that you cite the above article when using these data or modified code to prepare a publication.`

The files contained by this repository are numbered sequentially in the order they appear in our data analysis and figure generation workflow, each of which is described below. To use our code, you will need R installed with the AICcmodavg, broom, ggplot2, nls2, and reshape2 libraries, including their dependencies. You will also need to install the multifunc library from Jarrett Byrnes (http://github.com/jebyrnes/multifunc).


## `1 mean plot-level responses.csv`

The plot-level data used to generate all analyses and figures presented in the paper. These data represent the time-averaged value of each variable measured monthly in each plot over the course of the 10 month experiment. See the paper for a detailed description of the experiment and methodologies used. A full description of the methodologies and variables is also available from http://www.bco-dmo.org/dataset/716208. A brief description of each variable is given below. The suffix `.sr` denotes supporting responses only measured near the end of the experiment. 

`Plot` the experimental plot, numbered from West to East

`TrtPeg` the density treatment in total number of u-pegs assigned to each plot

`Gcvr` the average *Gracilaria* cover (%) maintained in each plot over the course of the experiment

`Epi` mean abundance of epifauna (# m^-2)

`EpiRich` mean richness of epifauna taxa (taxa m^-2)

`Dsln` mean chalk dissolution expressed as mass lost in grams per day (g d^-1)

`DslnFlip` mean reflected chalk dissolution

`Sed` mean sediment stabilization expressed as the change in height in cm per month (∆cm mo^-1)

`Nrsy` mean abundance of nursery species (# m^-2)

`NrsyRich` mean richness of nursery taxa (taxa m^-2)

`Dcmp` mean decomposition of *Spartina* stems expressed as mass lost per month (g mo^-1)

`Infa.sr` mean abundance of infauna (# L^-1)

`InfaRich.sr` mean richness of infauna taxa (taxa L^-1)

`Rays.sr` mean number of ray holes (# m^-2 d^-1)

`RaysFlip.sr` reflected mean number of ray holes

`Wfwl.sr` mean abundance of waterfowl (# m^-2 h^-1)


## `2 supplementary code.R `

Code to replicate the primary analysis of ecosystem multifunctionality as presented in Ramus et al. (2017) PNAS. This code calculates multifunctionaity, fits candidate models using nonlinear least squares for each function individually, and generates the model selection table and corresponding figure presented in the paper.


## `3 nlsTab.R `

Code to generate a model selection table for fitted nls models. This code works as a 'manual loop' for lack of a better description. Much of it is dedicated to combining and reorganizing information, although there are a few calculations.


## `4 model selection table.csv`

Model selection table of fitted nls models generated by the code.


## `5 figure.png`

Figure generated by the code. Corresponds to Fig. 2 presented in the paper.

![5%20figure](5%20figure.png)
