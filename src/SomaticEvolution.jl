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


export Cell,
MultiSimulation,
SimulationTracker,
SimulationInput,
MultilevelInput,
BranchingMoranInput,
BranchingInput,
MoranInput,
MultilevelInput,
Simulation,
SampledData,
PlotVAF,
PlotInverseVAF,
VAFResult,

multilevel_simulation,
multilevel_simulation_fast,
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
run_multilevel_from_file


include("definetypes.jl")
include("simulations.jl")
include("multilevel.jl")
include("process.jl")
include("sampling.jl")
include("analyse.jl")
include("multirun.jl")
include("plot.jl")
include("util.jl")
include("statistics.jl")
include("io.jl")
include("trees.jl")



end