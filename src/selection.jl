"""
    AbstractSelection

Type that defines how and when fit mutants will arise in the simulation.

# Subtypes
- `NeutralSelection`
- `SelectionPredefined`
- `SelectionDistribution`
"""
abstract type AbstractSelection end

"""
    NeutralSelection <: AbstractSelection

Neutral selection type: all mutations have the same fitness.
"""
struct NeutralSelection <: AbstractSelection
end

"""
    SelectionPredefined <: AbstractSelection

Non-neutral selection type: fit mutations arise at given times and with given selection 
coefficients.

# Fields
- `mutant_selection::Vector{Float64}` -- selection coefficient for each fit mutant that 
    arises
- `mutant_selection::Vector{Int64}` -- time that each fit mutant arises
"""
struct SelectionPredefined <: AbstractSelection
    mutant_selection::Vector{Float64}
    mutant_time::Vector{Float64}
end

"""
    SelectionDistribution{D<:Distribution} <: AbstractSelection

Non-neutral selection type: fit mutations arise at each cell division with a given 
    probability and with selection coefficient drawn from the distribution, up to a maximum 
    number.

# Fields
- `distribution::D` -- used to draw selection coefficients using `rand([rng, ]distibution)`
- `mutation_probability::Float64` -- probability that a fit mutant arises during cell 
    division
- `maximum_subclones::Int64` -- maximum number of fit mutants that can occur

"""
struct SelectionDistribution{D<:Distribution} <:AbstractSelection
    distribution::D
    mutation_probability::Float64
    maximum_subclones::Int64
end

"""
    getmaxsubclones(selection::AbstractSelection)

Returns the maximum number of subclones that can be generated from `selection`.
"""
function getmaxsubclones end

getmaxsubclones(::NeutralSelection) = 1
getmaxsubclones(selection::SelectionPredefined) = length(selection.mutant_selection) + 1
getmaxsubclones(selection::SelectionDistribution) = selection.maximum_subclones

"""
    newsubclone_ready(selection::AbstractSelection, nsubclonescurrent, nsubclones, t, rng)

Returns `true` if a new fit mutant (subclone) should be introduced.
"""
function newsubclone_ready end

newsubclone_ready(::NeutralSelection, nsubclonescurrent, nsubclones, t, rng) = false

function newsubclone_ready(selection::SelectionPredefined, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones && t >= selection.mutant_time[nsubclonescurrent]
end

function newsubclone_ready(selection::SelectionDistribution, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones &&  rand(rng) < selection.mutation_probability 
end

"""
    getselectioncoefficient(selection::AbstractSelection, nsubclonescurrent, rng)
 
Returns the next selection coefficient for a new mutant.
"""
function getselectioncoefficient end

function getselectioncoefficient(selection::SelectionPredefined, nsubclonescurrent, rng)
    selection.mutant_selection[nsubclonescurrent]
end
function getselectioncoefficient(selection::SelectionDistribution, nsubclonescurrent, rng)
    return rand(rng, selection.distribution)
end



"""
    Define a new selection type with all necessary functions

    Generate a new subclone writh mutation probability during division 
with new birth and death rate defined here. All clones are behaving in the same way.
"""
struct SelectionUniform <: AbstractSelection
    mutant_birthrate::Float64
    mutant_deathrate::Float64
    mutation_probability::Float64
    maximum_subclones::Int64
end

getmaxsubclones(selection::SelectionUniform) = selection.maximum_subclones

function newsubclone_ready(selection::SelectionUniform, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones &&  rand(rng) < selection.mutation_probability 
end
