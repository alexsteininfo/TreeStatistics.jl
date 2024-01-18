struct SampledData
    VAF::Vector{Float64}
    depth::Vector{Int64}
end

abstract type SimulationResult end

struct Simulation{S<:SimulationInput, T<:AbstractPopulation} <: SimulationResult
    input::S
    output::T
end

const MultiSimulation = Simulation{S, T} where {S<:MultilevelInput, T<:Population}

Base.length(multisim::MultiSimulation) = length(multisim.output)
Base.iterate(multisim::MultiSimulation) = iterate(multisim.output)
Base.iterate(multisim::MultiSimulation, state) = iterate(multisim.output, state)
Base.getindex(multisim::MultiSimulation, i) = getindex(multisim.output, i)
Base.setindex(multisim::MultiSimulation, v, i) = getindex(multisim.output, v, i)
Base.firstindex(multisim::MultiSimulation) = firstindex(multisim.output)
Base.lastindex(multisim::MultiSimulation) = lastindex(multisim.output)

get_simulation(multsim, i) = return Simulation(multsim.input, multsim.output[i])

abstract type AbstractVAFResult end
struct VAFResult{T<:SimulationInput} <: AbstractVAFResult
    read_depth::Float64
    cellularity::Float64
    detectionlimit::Float64
    trueVAF::Vector{Float64}
    sampledVAF::Vector{Float64}
    subclonefreq::Vector{Float64}
    subclonefreqp::Vector{Float64}
end

struct VAFResultMulti{T<:MultilevelInput} <: AbstractVAFResult
    read_depth::Float64
    cellularity::Float64
    detectionlimit::Float64
    trueVAFs::Vector{Vector{Float64}}
    sampledVAFs::Vector{Vector{Float64}}
    subclonefreqs::Vector{Vector{Float64}}
    subclonefreqps::Vector{Vector{Float64}}

end

"""
    show(io::IO, simulation::Simulation)

Print out summary of simulation.
"""
function Base.show(io::IO, simulation::SimulationResult)
    @printf(io, "===================================================================\n")
    show(io, simulation.input)
    # @printf(io, "\n")
    @printf(io, "===================================================================\n")
    show(io, simulation.output)
    @printf(io, "\n===================================================================")
end