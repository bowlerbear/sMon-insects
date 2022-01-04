# sMon - insects

A project to analyze heterogeneous data in Germany to examine the trends of dragonflies over the last decades

The repository contains the following folders:

derived-data: folder of the processed dragonfly data ready for analysis (including the archived folder of the data used for the Diversity and Distributions paper)

environ-data: scripts to process different data associated with the TK25 grid including naturraume, protected areasm land use etc..

formatting: formatting file for each dragonfly dataset obtained from each data provider (usuaully each federal state)

HPC-scripts: R scripts for the main analyses pushed to the HPC cluster - all the main analysis is here
For Diversity and Distributions opaper, the script is analysis_HPC_nation_naturruaum_sparta.R

model-auxfiles: associated data and helper files used by the models, e.g, task id text files linking tasks to jobs

old: bunch of files that were used but became redudant at some point. Best to ignore.

plots: plots for the papers and of the spline maps

R: R scripts of analysis run on the local PC e.g, exmaining the outputs of models run on the HPC, as well as some helper function files

raw-data: data files received by the data providers

splines: code to fit spline models. Note the models are processed in a separate repository - distributionChance - since they will form the next paper

state-analysis: analysis for specific states that were requested

trait-data: formatting the trait files. Note: the raw trait files are kept outside of the repository.
