abstract type AbstractCell end

"""
    Cell 

Represents a single cell.
"""
mutable struct Cell <: AbstractCell
    mutations::Vector{Int64}
    clonetype::Int64
    birthtime::Float64
    id::Int64
    parentid::Int64
end

Cell(mutations, clonetype) = Cell(mutations, clonetype, 0.0, 0, 0)
Cell(mutations, clonetype, birthtime) = Cell(mutations, clonetype, birthtime, 0, 0)

abstract type AbstractTreeCell <: AbstractCell end
mutable struct TreeCell <: AbstractTreeCell
    id::Int64
    alive::Bool
    birthtime::Float64
    mutations::Int64
    clonetype::Int64
end

TreeCell(; id=1, alive=true, birthtime=0.0, mutations=0, clonetype=1) =
    TreeCell(id, alive, birthtime, mutations, clonetype)

mutable struct SimpleTreeCell <: AbstractTreeCell
    id::Int64
    birthtime::Float64
    mutations::Int64
    clonetype::Int64
end

SimpleTreeCell(; id=1, birthtime=0.0, mutations=0, clonetype=1) =
    SimpleTreeCell(id, birthtime, mutations, clonetype)

function Base.show(io::IO, cell::TreeCell)
    print(io, "($(cell.id)) mutations = $(cell.mutations), t = $(cell.birthtime)")
    cell.alive || print(io, " X")
end

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

abstract type AbstractModule end
struct CellModule <: AbstractModule
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


Base.length(abstractmodule::AbstractModule) = abstractmodule.Nvec[end]
