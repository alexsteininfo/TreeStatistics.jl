struct SampledData
    VAF::Vector{Float64}
    depth::Vector{Int64}
end

abstract type SimulationResult end

struct Simulation{S<:SimulationInput, T<:AbstractModule} <: SimulationResult
    input::S
    output::T
end

struct MultiSimulation{S<:SimulationInput, T<:AbstractModule} <: SimulationResult
    input::S
    output::Vector{T}
end

Base.length(multisim::MultiSimulation) = length(multisim.output)
Base.iterate(multisim::MultiSimulation) = iterate(multisim.output)
Base.iterate(multisim::MultiSimulation, state) = iterate(multisim.output, state)
Base.getindex(multisim::MultiSimulation, i) = getindex(multisim.output, i)
Base.setindex(multisim::MultiSimulation, v, i) = getindex(multisim.output, v, i)
Base.firstindex(multisim::MultiSimulation) = firstindex(multisim.output)
Base.lastindex(multisim::MultiSimulation) = lastindex(multisim.output)

get_simulation(multsim, i) = return Simulation(multsim.input, multsim.output[i])

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

struct VAFResultMulti{T<:MultilevelInput}
    read_depth::Float64
    cellularity::Float64
    detectionlimit::Float64
    input::T
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
    @printf(io, "\n")
    show(io, simulation.output)

end

function Base.show(io::IO, abstractmodule::AbstractModule)
    @printf(io, "Final size = %d cells\n", length(abstractmodule))
    @printf(io, "Final time = %.2f", age(abstractmodule))
end

function Base.show(io::IO, population::Vector{T}) where T<:AbstractModule
    if length(population) == 0
        @printf(io, "Empty population")
    else
        @printf(io, "%d-module Population{%s} at time %.2f.\n", length(population), T, age(population))
        @printf(io, "Module sizes:")
        show(io, [length(m) for m in population])
    end
end

function Base.show(io::IO, input::BranchingInput)
    @printf(io, "Single level branching process:\n")
    @printf(io, "    Maximum cells = %d\n", input.Nmax)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Birth rate = %.2f, death rate = %.2f (branching)\n", input.birthrate, input.deathrate)
    @printf(io, "    %s: μ = %.1f\n", mutationdist_string(input.mutationdist), input.μ)
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
end

function Base.show(io::IO, input::MoranInput)
    @printf(io, "Single level Moran process:\n")
    @printf(io, "    Maximum cells = %d\n", input.N)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Birth/death rate = %.2f (Moran)\n", input.moranrate)
    @printf(io, "    %s: μ = %.1f\n", mutationdist_string(input.mutationdist), input.μ)
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
end

function Base.show(io::IO, input::BranchingMoranInput)
    @printf(io, "Single level Branching -> Moran process:\n")
    @printf(io, "    Maximum cells = %d\n", input.Nmax)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Birth rate = %.2f, death rate = %.2f (branching)\n", input.birthrate, input.deathrate)
    @printf(io, "    Birth/death rate = %.2f (Moran)\n", input.moranrate)
    @printf(io, "    %s: μ = %.1f\n", mutationdist_string(input.mutationdist), input.μ)
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
end

function Base.show(io::IO, input::MultilevelBranchingInput)
    @printf(io, "Multilevel branching process:\n")
    @printf(io, "    Maximum modules = %d\n", input.maxmodules)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Branch rate = %.2f\n", input.branchrate)
    @printf(io, "    Module size = %d\n", input.modulesize)
    @printf(io, "    Initial module size = %d\n", input.branchinitsize)
    @printf(io, "    Birth rate = %.2f, death rate = %.2f (branching)\n", input.birthrate, input.deathrate)
    @printf(io, "    Birth/death rate = %.2f (Moran)\n", input.moranrate)
    @printf(io, "    %s: μ = %.1f\n", mutationdist_string(input.mutationdist), input.μ)
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
end

function Base.show(io::IO, input::MultilevelBranchingMoranInput)
    @printf(io, "Multilevel branching -> Moran process:\n")
    @printf(io, "    Maximum modules = %d\n", input.maxmodules)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Branch rate = %.2f\n", input.branchrate)
    @printf(io, "    Module size = %d\n", input.modulesize)
    @printf(io, "    Initial module size = %d\n", input.branchinitsize)
    @printf(io, "    Birth rate = %.2f, death rate = %.2f (branching)\n", input.birthrate, input.deathrate)
    @printf(io, "    Birth/death rate = %.2f (Moran)\n", input.moranrate)
    @printf(io, "    %s: μ = %.1f\n", mutationdist_string(input.mutationdist), input.μ)
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
end

function mutationdist_string(mutationdist)
    if mutationdist == :fixed
        return "Fixed mutations"
    elseif mutationdist == :fixedtimedep
        return "Time-dependent fixed mutations"
    elseif mutationdist == :poisson
        return "Poisson distributed mutations"
    elseif mutationdist == :poissontimedep
        return "Time-dependent Poisson distibuted mutations"
    elseif mutationdist == :geometric
        return "Geometric distributed mutations"
    end
end