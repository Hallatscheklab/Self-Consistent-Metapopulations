Self-Consistent-Metapopulations

Numeric simulation of a Metacommunity undergoing Lotka-Volterra interactions with migration and demographic noise

Example code and evaluation used in this paper: https://www.pnas.org/doi/10.1073/pnas.2200390119


MetacommunityDynamics.py: includes an example Python code to solve the stochastic metacommunity dynamics as detailed in the SI Appendix, section 1 (we used Python version 3.6.3). The parameters, including number of species etc., can be set when running the code in python, e.g. type
"python MetacommunityDynamics.py TotalTime timestep Seed P alpha r S LogOfLambda K". The code generates a csv-file, which can be evaluated for instance using Mathematica (we used Version 12). The type of migration (short-range or global) can be chosen in the funtion update. Differences in the migration, interaction and growth rates among the species can be included by choosing rstd, alphstd, and diffstd nonzero in the code. Running the code with default parameters (i.e. "python MetacommunityDynamics") yields example results for a system of 100 species where species are close to their extinction threshold similar as shown in Fig.3 of the main text of the paper. 

ExampleEvaluation.nb: Mathematica notebook, which includes example code to plot the species' abundance profiles over time, the mean abundances of all demes over time, the mean abundance of all species over time, and the species' kymographs. In analogy, all figures in the main text and the Appendix SI can be obtained by evaluating the data as explained in the paper. 



