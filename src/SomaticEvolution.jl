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


export Cell,
MultiSimulation,
SimulationTracker,
SimulationInput,
MultilevelInput,
BranchingMoranInput,
BranchingInput,
MoranInput,
MultilevelInput,
InputParameters,
Simulation,
SampledData,
PlotVAF,
PlotInverseVAF,

multilevel_simulation,
multiplesimulations,
run1simulation, 
sampledhist,
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
pairwise_fixed_differences

include("definetypes.jl")
include("simulations.jl")
include("multilevel.jl")
include("process.jl")
include("sampling.jl")
include("analyse.jl")
include("multirun.jl")
include("plot.jl")
include("util.jl")
include("population_analysis.jl")



end