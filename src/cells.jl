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