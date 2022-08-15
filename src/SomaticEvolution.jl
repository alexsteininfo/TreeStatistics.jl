module SomaticEvolution

using Distributions
using Statistics
using Random
using StatsBase
using DataFrames
using RecipesBase
using GLM
using LaTeXStrings
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
SimulationTracker,
ModuleTracker,
TreeModuleTracker,
SimulationInput,
MultilevelInput,
BranchingMoranInput,
BranchingInput,
MoranInput,
MultilevelInput,
MultilevelBranchingInput,
MultilevelMoranInput,
Simulation,
SampledData,
PlotVAF,
PlotInverseVAF,
VAFResult,

multilevel_simulation,
multilevel_simulation_fast,
multilevel_simulation_timeseries,
multiplesimulations,
run1simulation, 
getallelefreq,
sampledallelefreq,
fitinverse,
gethist,
plotvaf,
plotvaf!,
plotinversevaf,
plotinversevaf!,
getstats,
saveinput,
branchingprocess,
moranprocess,
get_simulation,
pairwisedistance,
pairwise_fixed_differences,
pairwise_fixed_differences_matrix,
pairwise_fixed_differences_statistics,
shared_fixed_mutations,
newmoduletimes,
numbermodules,
cellpopulationsize,
meanmodulesize,
average_mutations,
average_mutations_per_module,
mutations_per_cell,
mutation_ids_by_cell,
clonal_mutation_ids,
clonal_mutations,
subclonefreq,
getVAFresult,
age,
saveinput,
loadinput,
run_multilevel_from_file,
endtime,
celllifetime,
celllifetimes,
getalivecells,
leftchild!,
rightchild!,
initialize,
run1simulation_tree,
time_to_MRCA,
coalescence_times,
getsingleroot,
popsize,
newinput


include("trees.jl")
include("cells.jl")
include("definetypes.jl")
include("simulations.jl")
include("multilevel.jl")
include("multilevel_trees.jl")
include("process.jl")
include("sampling.jl")
include("analyse.jl")
include("multirun.jl")
include("plot.jl")
include("util.jl")
include("statistics.jl")
include("io.jl")
include("simulation_trees.jl")




end