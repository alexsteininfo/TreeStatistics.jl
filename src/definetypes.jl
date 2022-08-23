
mutable struct CloneTracker
    parenttype::Int64
    parentmodule::Int64
    time::Float64
    mutations::Vector{Int64}
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

abstract type AbstractModule end
struct ModuleTracker <: AbstractModule
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
    cells::Vector{Cell}
    subclones::Vector{CloneTracker}
    id::Int64
    parentid::Int64
end

struct TreeModule{T<:AbstractTreeCell} <: AbstractModule 
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
    cells::Vector{BinaryNode{T}}
    subclones::Vector{CloneTracker}
    id::Int64
    parentid::Int64
end


Base.length(moduletracker::AbstractModule) = moduletracker.Nvec[end]

abstract type SimulationInput end

struct BranchingInput <: SimulationInput
    numclones::Int64 
    Nmax::Int64
    tmax::Float64
    clonalmutations::Int64
    selection::Vector{Float64}
    μ::Float64
    b::Float64
    d::Float64
    tevent::Vector{Float64}
    mutationdist::Symbol
    maxclonesize::Union{Int64, Nothing}
    ploidy::Int64
end

struct MoranInput <: SimulationInput
    N::Int64
    numclones::Int64 
    tmax::Float64
    clonalmutations::Int64
    selection::Vector{Float64}
    μ::Float64
    bdrate::Float64
    tevent::Vector{Float64}
    mutationdist::Symbol
    ploidy::Int64
end

struct BranchingMoranInput <: SimulationInput
    numclones::Int64 
    Nmax::Int64
    tmax::Float64
    clonalmutations::Int64
    selection::Vector{Float64}
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    tevent::Vector{Float64}
    mutationdist::Symbol
    ploidy::Int64
end

abstract type MultilevelInput <: SimulationInput end
struct MultilevelBranchingInput <: MultilevelInput
    modulesize::Int64
    maxtime::Float64
    maxmodules::Int64
    clonalmutations::Int64
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    mutationdist::Symbol
    branchrate::Float64
    branchinitsize::Int64
    ploidy::Int64
end

struct MultilevelMoranInput <: MultilevelInput
    modulesize::Int64
    maxtime::Float64
    maxmodules::Int64
    clonalmutations::Int64
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    mutationdist::Symbol
    branchrate::Float64
    branchinitsize::Int64
    ploidy::Int64
end
struct SimulationResult
    subclones::Vector{Clone}
    tend::Float64
    trueVAF::Vector{Float64}
    cells::Vector{Cell}
    Nvec::Vector{Int64}
    tvec::Vector{Float64}
end

struct SampledData
    VAF::Vector{Float64}
    depth::Vector{Int64}
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
    trueVAF::Vector{Float64}
    sampledVAF::Vector{Float64}
    subclonefreq::Vector{Float64}
    subclonefreqp::Vector{Float64}

end

struct MultiSimulation{S<:SimulationInput, T<:AbstractModule}
    input::S
    output::Vector{T}
end

# function MultiSimulation(::Type{T}, input::SimulationInput) where T <:AbstractModule
#     return MultiSimulation(
#         input,
#         T[]
#     )
# end

# MultiSimulation(input::SimulationInput) = MultiSimulation(input, ModuleTracker[])

Base.length(multisim::MultiSimulation) = length(multisim.output)
Base.iterate(multisim::MultiSimulation) = iterate(multisim.output)
Base.iterate(multisim::MultiSimulation, state) = iterate(multisim.output, state)
Base.getindex(multisim::MultiSimulation, i) = getindex(multisim.output, i)
Base.setindex(multisim::MultiSimulation, v, i) = getindex(multisim.output, v, i)
Base.firstindex(multisim::MultiSimulation) = firstindex(multisim.output)
Base.lastindex(multisim::MultiSimulation) = lastindex(multisim.output)



# const MultiSimulation  = MultiSimulation{T, ModuleTracker} where T<:MultilevelInput
# const MultiSimulation = MultiSimulation{T, TreeModule} where T<:MultilevelInput 


function get_simulation(multsim, i)
    return Simulation(multsim.input, multsim.output[i])
end 

function BranchingInput(;numclones=1, Nmax=10000, tmax=Inf, ploidy=2, μ=10.0, 
    clonalmutations=0, selection=fill(0.0,numclones), b=log(2.0), d=0.0, 
    tevent=collect(1.0:0.5:(1+numclones)/2), fixedmu=false, mutationdist=nothing,
    maxclonesize=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

    return BranchingInput(
        numclones,
        Nmax,
        tmax,
        clonalmutations,
        selection,
        μ,
        b,
        d,
        tevent,
        mutationdist,
        maxclonesize,
        ploidy
    )
end

function MoranInput(;numclones=1, N=10000, ploidy=2, μ=10.0, clonalmutations=0, 
    selection=fill(0.0,numclones), bdrate=log(2.0), tmax=15.0,
    tevent=collect(1.0:0.5:(1+numclones)/2), fixedmu=false, mutationdist=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

    return MoranInput(
        N,
        numclones,
        tmax,
        clonalmutations,
        selection,
        μ,
        bdrate,
        tevent,
        mutationdist,
        ploidy
    )
end

function BranchingMoranInput(;numclones=1, Nmax=10000, ploidy=2, μ=10.0, 
    clonalmutations=0, selection=fill(0.0,numclones), bdrate=log(2.0), b=log(2), 
    d=0, tmax=15.0, tevent=collect(1.0:0.5:(1+numclones)/2), fixedmu=false, 
    mutationdist=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

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
        mutationdist,
        ploidy
    )
    
end

function MultilevelBranchingInput(;modulesize=200, ploidy=2, μ=10.0, clonalmutations=0, 
    bdrate=1, b=1, d=0, maxtime=15, maxmodules=10000, fixedmu=false, 
    mutationdist=nothing, branchrate=5, branchfraction=0.1, branchinitsize=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

    return MultilevelBranchingInput(
            modulesize,
            maxtime,
            maxmodules,
            clonalmutations,
            μ,
            bdrate,
            b,
            d,
            mutationdist,
            branchrate,
            branchinitsize !== nothing ? branchinitsize : ceil(modulesize * branchfraction),
            ploidy,
    )
end

function MultilevelMoranInput(;modulesize=200, ploidy=2, μ=10.0, clonalmutations=0, 
    bdrate=1, b=1, d=0, maxtime=15, maxmodules=10000, fixedmu=false, 
    mutationdist=nothing, branchrate=0.1, branchfraction=0.1, branchinitsize=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

    return MultilevelMoranInput(
            modulesize,
            maxtime,
            maxmodules,
            clonalmutations,
            μ,
            bdrate,
            b,
            d,
            mutationdist,
            branchrate,
            branchinitsize !== nothing ? branchinitsize : ceil(modulesize * branchfraction),
            ploidy,
    )
end

function newinput(input::T; kwargs...) where T <: SimulationInput
    newkwargs = Dict(
        field in keys(kwargs) ? field => kwargs[field] : field => getfield(input, field)
            for field in fieldnames(T))
    return T(;kwargs...)
end

age(moduletracker::AbstractModule) = moduletracker.tvec[end]
age(populationtracker::Vector{T}) where T<:AbstractModule = maximum(map(age, populationtracker))
age(simulation::Simulation) = age(simulation.output)
age(multisim::MultiSimulation) = maximum(age(output) for output in multisim.output)

function set_mutationdist(mutationdist, fixedmu)
    if isnothing(mutationdist)
        return fixedmu ? :fixed : :poisson
    elseif typeof(mutationdist) === String
        return Symbol(mutationdist)
    else 
        return mutationdist
    end
end