"""
    Cell 

Represents a single cell.
"""
mutable struct Cell
    mutations::Array{Int64,1}
    clonetype::Int64
end

struct CloneTracker
    parenttype::Int64
    time::Float64
    mutations::Array{Int64, 1}
    N0::Int64
    Ndivisions::Int64
    avdivisions::Float64
end

struct Clone
    parenttype::Int64
    time::Float64
    mutations::Int64
    N0::Int64
    Ndivisions::Int64
    avdivisions::Float64
    freq::Float64
    freqp::Float64
end

abstract type SimulationTracker end

struct BranchingTracker <: SimulationTracker
    Nvec::Array{Int64, 1}
    tvec::Array{Float64, 1}
    cells::Array{Cell, 1}
    birthrates::Array{Float64, 1}
    deathrates::Array{Float64, 1}
    clonesize::Array{Int64, 1}
    subclones::Array{CloneTracker, 1}
end

struct MoranTracker <: SimulationTracker
    N::Int64
    tvec::Array{Float64, 1}
    cells::Array{Cell, 1}
    clonesize::Array{Int64, 1}
    subclones::Array{CloneTracker, 1}
end

abstract type SimulationInput end

struct BranchingInput <: SimulationInput
    numclones::Int64 
    Nmax::Int64
    clonalmutations::Int64
    selection::Array{Float64,1}
    μ::Float64
    b::Float64
    d::Float64
    tevent::Array{Float64,1}
    fixedmu::Bool
    maxclonesize::Int64
end

struct MoranInput <: SimulationInput
    N::Int64
    numclones::Int64 
    tmax::Float64
    clonalmutations::Int64
    selection::Array{Float64,1}
    μ::Float64
    bdrate::Float64
    tevent::Array{Float64,1}
    fixedmu::Bool
end

struct InputParameters{T<:SimulationInput}
    detectionlimit::Float64
    ploidy::Int64
    read_depth::Float64
    ρ::Float64
    cellularity::Float64
    siminput::T
end

struct SimulationResult
    subclones::Array{Clone, 1}
    tend::Float64
    trueVAF::Array{Float64,1}
    cells::Array{Cell, 1}
end

struct SampledData
    VAF::Array{Float64,1}
    depth::Array{Int64,1}
end

struct Simulation{T<:SimulationInput}
    input::InputParameters{T}
    output::SimulationResult
    sampled::SampledData
end

struct MultiSimulation{T<:SimulationInput}
    input::InputParameters{T}
    output::Array{SimulationResult,1}
    sampled::Array{SampledData,1}
end

function get_simulation(multsim::MultiSimulation, i)
    return Simulation(multsim.input, multsim.output[i], multsim.sampled[i])
end 

function InputParameters{BranchingInput}(;numclones = 1, Nmax = 10000, ploidy = 2, 
    read_depth = 100.0, detectionlimit = 5/read_depth, μ = 10.0, clonalmutations = μ, 
    selection = fill(0.0,numclones), b = log(2.0), d = 0.0, 
    tevent = collect(1.0:0.5:(1+numclones)/2), ρ = 0.0, cellularity = 1.0, fixedmu = false, 
    maxclonesize = 200)

    return InputParameters(
        detectionlimit,
        ploidy,
        read_depth,
        ρ,
        cellularity,
        BranchingInput(
            numclones,
            Nmax,
            clonalmutations,
            selection,
            μ,
            b,
            d,
            tevent,
            fixedmu,
            maxclonesize
        )
    )
end

function InputParameters{MoranInput}(;numclones = 1, N = 10000, ploidy = 2, 
    read_depth = 100.0, detectionlimit = 5/read_depth, μ = 10.0, clonalmutations = μ, 
    selection = fill(0.0,numclones), bdrate = log(2.0), tmax = 15.0,
    tevent = collect(1.0:0.5:(1+numclones)/2), ρ = 0.0, cellularity = 1.0, fixedmu = false)

    return InputParameters(
        detectionlimit,
        ploidy,
        read_depth,
        ρ,
        cellularity,
        MoranInput(
            N,
            numclones,
            tmax,
            clonalmutations,
            selection,
            μ,
            bdrate,
            tevent,
            fixedmu,
        )
    )
end
