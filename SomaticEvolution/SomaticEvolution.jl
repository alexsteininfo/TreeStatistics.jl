module SomaticEvolution

using Distributions
using Statistics
using Distributions
using Random
using StatsBase
using DataFrames

export Cell,
SimulationTracker,
SimulationInput,
BranchingInput,
MoranInput,
InputParameters,
Simulation,
SampledData,

run1simulation, 
sampledhist,
gethist

include("runsimulations.jl")
include("process.jl")
include("sampling.jl")



end