#region Define BinaryNode and methods for binary trees

"""
    BinaryNode{T}

Basic unit of a binary tree.

# Fields
- `data::T`
- `parent::Union{Nothing, BinaryNode{T}}`
- `left::Union{Nothing, BinaryNode{T}}`
- `right::Union{Nothing, BinaryNode{T}}`
"""
mutable struct BinaryNode{T}
    data::T
    parent::Union{Nothing, BinaryNode{T}}
    left::Union{Nothing, BinaryNode{T}}
    right::Union{Nothing, BinaryNode{T}}

    function BinaryNode{T}(data, parent=nothing, l=nothing, r=nothing) where T
        new{T}(data, parent, l, r)
    end
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

"""
    leftchild!(parent::BinaryNode, data)

Create a new `BinaryNode` from `data` and assign it to `parent.left`.

See also [`rightchild!`](@ref).
"""
function leftchild!(parent::BinaryNode, data)
    isnothing(parent.left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

"""
    rightchild!(parent::BinaryNode, data)

Create a new `BinaryNode` from `data` and assign it to `parent.right`.

See also [`leftchild!`](@ref).
"""
function rightchild!(parent::BinaryNode, data)
    isnothing(parent.right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end

function AbstractTrees.children(node::BinaryNode)
    if isnothing(node.left) && isnothing(node.right)
        ()
    elseif isnothing(node.left) && !isnothing(node.right)
        (node.right,)
    elseif !isnothing(node.left) && isnothing(node.right)
        (node.left,)
    else
        (node.left, node.right)
    end
end

function AbstractTrees.nextsibling(child::BinaryNode)
    isnothing(child.parent) && return nothing
    p = child.parent
    if !isnothing(p.right)
        child === p.right && return nothing
        return p.right
    end
    return nothing
end

function AbstractTrees.prevsibling(child::BinaryNode)
    isnothing(child.parent) && return nothing
    p = child.parent
    if !isnothing(p.left)
        child === p.left && return nothing
        return p.left
    end
    return nothing
end

AbstractTrees.nodevalue(n::BinaryNode) = n.data
AbstractTrees.ParentLinks(::Type{<:BinaryNode}) = StoredParents()
AbstractTrees.parent(n::BinaryNode) = n.parent
AbstractTrees.NodeType(::Type{<:BinaryNode{T}}) where {T} = HasNodeType()
AbstractTrees.nodetype(::Type{<:BinaryNode{T}}) where {T} = BinaryNode{T}

#ensure all nodes in tree are of the same type
Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()

AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

"""
    popsize(root::BinaryNode)

Get the number of alive cells that are descendents of `root`.
"""
function popsize(root::BinaryNode)
    N = 0
    for l in Leaves(root)
        if isalive(nodevalue(l))
            N += 1
        end
    end
    return N
end

haschildren(node::BinaryNode) = length(children(node)) != 0

function Base.show(io::IO, node::BinaryNode)
    show(io, node.data)
end

function Base.show(io::IO, nodevec::Vector{BinaryNode{T}}) where T
    println(io, "$(length(nodevec))-element Vector{BinaryNode{$T}}:")
    for node in nodevec
        show(io, node)
        print(io, "\n")
    end
end
#endregion

#region Define cell types
abstract type AbstractCell end

"""
    Cell

Represents a single cell with a list of mutation ids.
"""
mutable struct Cell <: AbstractCell
    mutations::Vector{Int64}
    clonetype::Int64
    birthtime::Float64
    latestupdatetime::Float64
    id::Int64
    parentid::Int64
end

Cell(mutations, clonetype) = Cell(mutations, clonetype, 0.0, 0.0, 0, 0)

function Cell(mutations, clonetype, birthtime)
    return Cell(mutations, clonetype, birthtime, birthtime, 0, 0)
end

abstract type AbstractTreeCell <: AbstractCell end
mutable struct TreeCell <: AbstractTreeCell
    id::Int64
    alive::Bool
    birthtime::Float64
    latestupdatetime::Float64
    mutations::Int64
    clonetype::Int64
end
"""
    TreeCell

Represents a single cell that can be the `data` field of a `BinaryNode{SimpleTreeCell}`.
Has an `alive` field, so dead cells are kept in tree structure.
[`SimpleTreeCell`](@ref) is preferred.
"""
function TreeCell(; id=1, alive=true, birthtime=0.0, mutations=0, clonetype=1)
    return TreeCell(id, alive, birthtime, birthtime, mutations, clonetype)
end

"""
    SimpleTreeCell

Represents a single cell that can be the `data` field of a `BinaryNode{SimpleTreeCell}`.
Does not have an `alive` field, so dead cells must be removed from tree
structure. Preferred over [`TreeCell`](@ref).
"""
mutable struct SimpleTreeCell <: AbstractTreeCell
    id::Int64
    birthtime::Float64
    latestupdatetime::Float64
    mutations::Int64
    clonetype::Int64
end

function SimpleTreeCell(; id=1, birthtime=0.0, mutations=0, clonetype=1)
    return SimpleTreeCell(id, birthtime, birthtime, mutations, clonetype)
end

function Base.show(io::IO, cell::TreeCell)
    print(io, "($(cell.id)) mutations = $(cell.mutations), t = $(cell.birthtime)")
    cell.alive || print(io, " X")
end

getclonetype(cell::Cell) = cell.clonetype
getclonetype(cellnode::BinaryNode) = cellnode.data.clonetype
setclonetype(cell::Cell, newclonetype) = (cell.clonetype = newclonetype)
setclonetype(cellnode::BinaryNode, newclonetype) = (cellnode.data.clonetype = newclonetype)
#endregion

"""
    Subclone

Defines subclone properties, including the `birthrate`, `deathrate`, `moranrate` and
`asymmetricrate` of cells within the subclone.
"""
@kwdef mutable struct Subclone
    subcloneid::Int64 = 1
    parentid::Int64 = 0
    mutationtime::Float64 = 0.0
    size::Int64 = 1
    birthrate::Float64 = 1.0
    deathrate::Float64 = 0.0
    moranrate::Float64 = 1.0
    asymmetricrate::Float64 = 0.0
end

Base.length(subclone::Subclone) = subclone.size

#region Define module structure types
"""
    ModuleStructure

Abstract type that can be used to implement different spatial structures
for modules. Currently only [`WellMixed`](@ref) is implemented.
"""
abstract type ModuleStructure end

"""
    WellMixed

ModuleStructure subtype that assumes no spatial structure.
"""
struct WellMixed <: ModuleStructure end

"""
    Linear

ModuleStructure subtype that assumes cells form a line within the module.
Not properly implemented.
"""
struct Linear <: ModuleStructure
    size::Int64
end
#endregion

#region Define module types and methods

"""
    AbstractModule

Abstract type that defines the module. Subtypes are [`CellModule`](@ref) and `TreeModule`(@ref).
"""
abstract type AbstractModule end

"""
    CellModule{S<:ModuleStructure} <: AbstractModule

Module type that stores cells as `Cell` type which each has a vector of mutation ids.
"""
mutable struct CellModule{S<:ModuleStructure} <: AbstractModule
    cells::Vector{Union{Cell, Nothing}}
    t::Float64
    branchtimes::Vector{Float64}
    id::Int64
    parentid::Int64
    structure::S
end

"""
    TreeModule{T<:AbstractTreeCell, S<:ModuleStructure} <: AbstractModule

Module type that stores cells as `SimpleTreeCell` or `TreeCell` type. Each cell the `data`
field of a `BinaryNode`, so the full ancestory is maintained as a tree structure.
"""
mutable struct TreeModule{T<:AbstractTreeCell, S<:ModuleStructure} <: AbstractModule
    cells::Vector{Union{BinaryNode{T}, Nothing}}
    t::Float64
    branchtimes::Vector{Float64}
    id::Int64
    parentid::Int64
    structure::S
end

const TreeCellVector = Union{
    Vector{Union{BinaryNode{TreeCell}, Nothing}},
    Vector{BinaryNode{TreeCell}}
}
const SimpleTreeCellVector = Union{
    Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}},
    Vector{BinaryNode{SimpleTreeCell}}
}
const AbstractTreeCellVector = Union{TreeCellVector, SimpleTreeCellVector}
const CellVector = Union{Vector{Union{Cell, Nothing}}, Vector{Cell}}
const AbstractCellVector = Union{CellVector, AbstractTreeCellVector}

"""
    moduletype(::Type{T}, ::Type{S})

Determine correct type for a new module where first type is a cell, vector of cells or binary nodes and
second type is a module structure.

"""
function moduletype end
moduletype(::Type{T}, ::Type{S}) where {T <: AbstractTreeCell, S} = TreeModule{T,S}
moduletype(::Type{BinaryNode{T}}, ::Type{S}) where {T <: AbstractTreeCell, S} = moduletype(T, S)
moduletype(::Type{Cell}, ::Type{S}) where S = CellModule{S}
moduletype(::Type{Vector{T}}, ::Type{S}) where {T, S} = moduletype(T, S)
moduletype(::Type{Vector{Union{T, Nothing}}}, ::Type{S}) where {T, S} = moduletype(T, S)

moduleid(abstractmodule::AbstractModule) = abstractmodule.id

function firstcellnode(treemodule::TreeModule)
    for cellnode in treemodule.cells
        if !isnothing(cellnode)
            return cellnode
        end
    end
end

allcells(abstractmodule) = filter(x -> !isnothing(x), abstractmodule.cells)

#endregion

#region Methods for finding most recent common ancestor and tree roots
"""
    getroot(nodevec)

Get a list of unique roots for each tree that contains a node in `nodevec`.
This requires node to have the trait StoredParents.
"""
AbstractTrees.getroot(nodevec::AbstractTreeCellVector) = collect(Set(getroot.(nodevec)))

"""
    getsingleroot(nodevec)

Find the unique root of a list of nodes `nodevec`. If there is no unique root return
`nothing`.
"""
function getsingleroot(nodevec::AbstractTreeCellVector)
    roots = AbstractTrees.getroot(nodevec)
    if length(roots) == 1
        return roots[1]
    else
        return nothing
    end
end

"""
    findMRCA(treemodule)
    findMRCA(nodes)
    findMRCA(node1, node2)

Find the most recent common ancestor of all nodes in `treemodule.cells`, list of nodes
`nodes` or two nodes `node1` and `node2`.
"""
function findMRCA end

function findMRCA(node1, node2)
    node1 == node2 && return node1
    while true
        #ensure that node1 is higher (or equal) level in tree than node2
        if  node1.data.id > node2.data.id
            node1, node2 = node2, node1
        end
        #if nodes are siblings return parent
        if node1.parent == node2.parent
            return node1.parent
        #if both nodes are roots there is no MRCA
        elseif isnothing(node1.parent) && isnothing(node2.parent)
            return nothing
        #if node1 is parent of node2 return node1
        elseif node1 == node2.parent
            return node1
        #if none of these are satisfied keep looking
        else
            node2 = node2.parent
        end
    end
end

function findMRCA(nodes::Vector)
    nodes = copy(nodes)
    node1 = pop!(nodes)
    while length(nodes) > 0
        node2 = pop!(nodes)
        node1 = findMRCA(node1, node2)
    end
    return node1
end


function findMRCA(treemodule)
    return findMRCA(treemodule.cells)
end
#endregion

Base.length(abstractmodule::AbstractModule) = length(abstractmodule.cells)
Base.iterate(abstractmodule::AbstractModule) = iterate(abstractmodule.cells)
Base.iterate(abstractmodule::AbstractModule, state) = iterate(abstractmodule.cells, state)
Base.getindex(abstractmodule::AbstractModule, i) = getindex(abstractmodule.cells, i)
Base.setindex(abstractmodule::AbstractModule, v, i) = getindex(abstractmodule.cells, v, i)
Base.firstindex(abstractmodule::AbstractModule) = firstindex(abstractmodule.cells)
Base.lastindex(abstractmodule::AbstractModule) = lastindex(abstractmodule.cells)
