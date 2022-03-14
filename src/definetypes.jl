
"""
    Cell 

Represents a single cell.
"""
mutable struct Cell
    mutations::Array{Int64,1}
    clonetype::Int64
end

mutable struct CloneTracker
    parenttype::Int64
    parentmodule::Int64
    time::Float64
    mutations::Array{Int64, 1}
    N0::Int64
    Ndivisions::Int64
    avdivisions::Float64
    size::Int64
end

struct Clone
    parenttype::Int64
    parentmodule::Int64
    time::Float64
    mutations::Int64
    N0::Int64
    Ndivisions::Int64
    avdivisions::Float64
    freq::Float64
    freqp::Float64
end
struct ModuleTracker 
    Nvec::Array{Int64, 1}
    tvec::Array{Float64, 1}
    cells::Array{Cell, 1}
    subclones::Array{CloneTracker, 1}
    id::Int64
    parentid::Int64
end

Base.length(moduletracker::ModuleTracker) = moduletracker.Nvec[end]

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
    ploidy::Int64
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
    ploidy::Int64
end

struct BranchingMoranInput <: SimulationInput
    numclones::Int64 
    Nmax::Int64
    tmax::Float64
    clonalmutations::Int64
    selection::Array{Float64,1}
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    tevent::Array{Float64,1}
    fixedmu::Bool
    ploidy::Int64
end
struct MultilevelInput <: SimulationInput
    modulesize::Int64
    maxtime::Float64
    maxmodules::Int64
    clonalmutations::Int64
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    fixedmu::Bool
    branchrate::Float64
    branchinitsize::Int64
    ploidy::Int64
end

struct SimulationResult
    subclones::Array{Clone, 1}
    tend::Float64
    trueVAF::Array{Float64,1}
    cells::Array{Cell, 1}
    Nvec::Array{Int64, 1}
    tvec::Array{Float64, 1}
end

struct SampledData
    VAF::Array{Float64,1}
    depth::Array{Int64,1}
end

struct Simulation{T<:SimulationInput}
    input::T
    output::ModuleTracker
end

struct VAFResult{T<:SimulationInput}
    read_depth::Float64
    cellularity::Float64
    detectionlimit::Float64
    input::T
    trueVAF::Array{Float64,1}
    sampledVAF::Array{Float64,1}
    subclonefreq::Array{Float64, 1}
    subclonefreqp::Array{Float64, 1}

end

struct MultiSimulation{T<:SimulationInput}
    input::T
    output::Array{ModuleTracker, 1}
    # mutations::Array{Int64, 1}
end

function MultiSimulation{T}(input::T) where T<:SimulationInput
    return MultiSimulation(
        input,
        ModuleTracker[]
    )
end

Base.length(multisim::MultiSimulation) = length(multisim.output)
Base.iterate(multisim::MultiSimulation) = iterate(multisim.output)
Base.iterate(multisim::MultiSimulation, state) = iterate(multisim.output, state)
Base.getindex(multisim::MultiSimulation, i) = getindex(multisim.output, i)
Base.setindex(multisim::MultiSimulation, v, i) = getindex(multisim.output, v, i)
Base.firstindex(multisim::MultiSimulation) = firstindex(multisim.output)
Base.lastindex(multisim::MultiSimulation) = lastindex(multisim.output)



const Population = MultiSimulation{MultilevelInput}

function get_simulation(multsim, i)
    return Simulation(multsim.input, multsim.output[i])
end 

function BranchingInput(;numclones = 1, Nmax = 10000, ploidy = 2, μ = 10.0, 
    clonalmutations = μ, selection = fill(0.0,numclones), b = log(2.0), d = 0.0, 
    tevent = collect(1.0:0.5:(1+numclones)/2), fixedmu = false, 
    maxclonesize = 200)

    return BranchingInput(
            numclones,
            Nmax,
            clonalmutations,
            selection,
            μ,
            b,
            d,
            tevent,
            fixedmu,
            maxclonesize,
            ploidy
    )
end

function MoranInput(;numclones = 1, N = 10000, ploidy = 2, μ = 10.0, clonalmutations = μ, 
    selection = fill(0.0,numclones), bdrate = log(2.0), tmax = 15.0,
    tevent = collect(1.0:0.5:(1+numclones)/2), fixedmu = false)

    return MoranInput(
            N,
            numclones,
            tmax,
            clonalmutations,
            selection,
            μ,
            bdrate,
            tevent,
            fixedmu,
            ploidy
    )
end

function BranchingMoranInput(;numclones = 1, Nmax = 10000, ploidy = 2, μ = 10.0, 
    clonalmutations = μ, selection = fill(0.0,numclones), bdrate = log(2.0), b = log(2), 
    d = 0, tmax = 15.0, tevent = collect(1.0:0.5:(1+numclones)/2), fixedmu = false)

    return BranchingMoranInput(
            numclones,
            Nmax,
            tmax,
            clonalmutations,
            selection,
            μ,
            bdrate,
            b,
            d,
            tevent,
            fixedmu,
            ploidy
    )
    
end

function MultilevelInput(;modulesize=200, ploidy=2, μ=10.0, clonalmutations=μ, 
    bdrate=log(2.0), b=log(2), d=0, maxtime=15, maxmodules=10000, fixedmu=false, branchrate=5, 
    branchfraction=0.1, branchinitsize=nothing)

    return MultilevelInput(
            modulesize,
            maxtime,
            maxmodules,
            clonalmutations,
            μ,
            bdrate,
            b,
            d,
            fixedmu,
            branchrate,
            branchinitsize !== nothing ? branchinitsize : ceil(modulesize * branchfraction),
            ploidy,
    )
end

age(moduletracker::ModuleTracker) = moduletracker.tvec[end]
age(populationtracker::Vector{ModuleTracker}) = maximum(map(age, populationtracker))
age(simulation::Simulation) = age(simulation.output)
age(multisim::MultiSimulation) = maximum(age(output) for output in multisim.output)