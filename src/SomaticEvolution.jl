module SomaticEvolution

using Distributions
using Statistics
using Random
using StatsBase
using DataFrames
using RecipesBase
using GLM
using Printf
using CSV
using Plots
using JSON
using DelimitedFiles
using AbstractTrees

export 
Cell,
TreeCell,
SimpleTreeCell,
BinaryNode,
MultiSimulation,
CellModule,
TreeCellModule,
SimulationInput,
BranchingMoranInput,
BranchingInput,
MoranInput,
MultilevelBranchingInput,
MultilevelBranchingMoranInput,
Simulation,
SampledData,
VAFResult,
VAFMultiResult,

#functions for running multilevel or single-level simulations
runsimulation,
runsimulation_condfixtime,
runsimulation_condfixtime_to_nfixed,
runsimulation_timeseries,
runsimulation_timeseries_returnfinalpop,
run1simulation, 
multilevel_simulation,
multilevel_simulation_timeseries,
run_multilevel_from_file,
initialize,

#functions for looking at VAF distributions (only implemented for Cell based simulations)
getallelefreq,
sampledallelefreq,
fitinverse,
gethist,
plotvaf,
plotvaf!,
plotinversevaf,
plotinversevaf!,
getVAFresult,
getVAFresultmulti,


#other statistics currently implemented only for Cell based simulations
shared_fixed_mutations,
subclonefreq,

#other statistics (work for Cell and (Simple)TreeCell simulations)
pairwisedistance,
pairwise_fixed_differences,
pairwise_fixed_differences_matrix,
pairwise_fixed_differences_statistics,
pairwise_fixed_differences_clonal,
newmoduletimes,
numbermodules,
cellpopulationsize,
meanmodulesize,
average_mutations,
average_mutations_per_module,
mutations_per_cell,
clonal_mutations,
age,

#functions for TreeModule
endtime,
celllifetime,
celllifetimes,
getalivecells,
leftchild!,
rightchild!,
time_to_MRCA,
coalescence_times,
getsingleroot,
popsize,
findMRCA,
moduleid,

#util
newinput,
saveinput,
get_simulation,
saveinput,
loadinput

include("input.jl")
include("trees.jl")
include("cells_modules.jl")
include("results.jl")
include("simulations.jl")
include("multilevel.jl")
include("multilevel_condfixtime.jl")
include("multilevel_trees.jl")
include("process.jl")
include("samplingVAF.jl")
include("analyseVAF.jl")
include("multirun.jl")
include("plot.jl")
include("util.jl")
include("statistics.jl")
include("io.jl")
include("simulation_trees.jl")




end