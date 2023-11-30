abstract type AbstractSelection end

struct NeutralSelection <: AbstractSelection
end

struct SelectionPredefined <: AbstractSelection
    mutant_selection::Vector{Float64}
    mutant_time::Vector{Float64}
end

struct SelectionDistribution{D<:Distribution} <:AbstractSelection
    distribution::D
    mutation_probability::Float64
    maximum_subclones::Int64
end

getmaxsubclones(::NeutralSelection) = 1
getmaxsubclones(selection::SelectionPredefined) = length(selection.mutant_selection) + 1
getmaxsubclones(selection::SelectionDistribution) = selection.maximum_subclones

newsubclone_ready(::NeutralSelection, nsubclonescurrent, nsubclones, t, rng) = false

function newsubclone_ready(selection::SelectionPredefined, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones && t >= selection.mutant_time[nsubclonescurrent]
end

function newsubclone_ready(selection::SelectionDistribution, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones &&  rand(rng) < selection.mutation_probability 
end

function getselectioncoefficient(selection::SelectionPredefined, nsubclonescurrent, rng)
    selection.mutant_selection[nsubclonescurrent]
end
function getselectioncoefficient(selection::SelectionDistribution, nsubclonescurrent, rng)
    return rand(rng, selection.distribution)
end
