#region Define BinaryNode and methods for binary trees
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

function leftchild!(parent::BinaryNode, data)
    isnothing(parent.left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end
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

Represents a single cell.
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
Cell(mutations, clonetype, birthtime) = Cell(mutations, clonetype, birthtime, birthtime, 0, 0)

abstract type AbstractTreeCell <: AbstractCell end
mutable struct TreeCell <: AbstractTreeCell
    id::Int64
    alive::Bool
    birthtime::Float64
    latestupdatetime::Float64
    mutations::Int64
    clonetype::Int64
end

TreeCell(; id=1, alive=true, birthtime=0.0, mutations=0, clonetype=1) =
    TreeCell(id, alive, birthtime, birthtime, mutations, clonetype)

mutable struct SimpleTreeCell <: AbstractTreeCell
    id::Int64
    birthtime::Float64
    latestupdatetime::Float64
    mutations::Int64
    clonetype::Int64
end

SimpleTreeCell(; id=1, birthtime=0.0, mutations=0, clonetype=1) =
    SimpleTreeCell(id, birthtime, birthtime, mutations, clonetype)

function Base.show(io::IO, cell::TreeCell)
    print(io, "($(cell.id)) mutations = $(cell.mutations), t = $(cell.birthtime)")
    cell.alive || print(io, " X")
end
#endregion

getclonetype(cell::Cell) = cell.clonetype
getclonetype(cellnode::BinaryNode) = cellnode.data.clonetype
setclonetype(cell::Cell, newclonetype) = (cell.clonetype = newclonetype)
setclonetype(cellnode::BinaryNode, newclonetype) = (cellnode.data.clonetype = newclonetype)
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
abstract type ModuleStructure end

struct WellMixed <: ModuleStructure end
struct Linear <: ModuleStructure 
    size::Int64
end
#endregion

#region Define module types and methods
abstract type AbstractModule end
mutable struct CellModule{S<:ModuleStructure} <: AbstractModule
    cells::Vector{Union{Cell, Nothing}}
    t::Float64
    branchtimes::Vector{Float64}
    id::Int64
    parentid::Int64
    structure::S
end

mutable struct TreeModule{T<:AbstractTreeCell, S<:ModuleStructure} <: AbstractModule 
    cells::Vector{Union{BinaryNode{T}, Nothing}}
    t::Float64
    branchtimes::Vector{Float64}
    id::Int64
    parentid::Int64
    structure::S
end

const TreeCellVector = Union{Vector{Union{BinaryNode{TreeCell}, Nothing}}, Vector{BinaryNode{TreeCell}}}
const SimpleTreeCellVector = Union{Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}, Vector{BinaryNode{SimpleTreeCell}}}
const AbstractTreeCellVector = Union{TreeCellVector, SimpleTreeCellVector}
const CellVector = Union{Vector{Union{Cell, Nothing}}, Vector{Cell}}
const AbstractCellVector = Union{CellVector, AbstractTreeCellVector}

#determine correct type for a new module based on cell type or vector type
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
AbstractTrees.getroot(nodevec::AbstractTreeCellVector) = collect(Set(getroot.(nodevec)))

function getsingleroot(nodevec::AbstractTreeCellVector)
    roots = AbstractTrees.getroot(nodevec)
    if length(roots) == 1
        return roots[1]
    else
        return nothing
    end
end

function findMRCA(cellnode1, cellnode2)
    cellnode1 == cellnode2 && return cellnode1
    while true
        #ensure that cellnode1 is higher (or equal) level in tree than cellnode2
        if  cellnode1.data.id > cellnode2.data.id
            cellnode1, cellnode2 = cellnode2, cellnode1
        end
        #if cell nodes are siblings return parent
        if cellnode1.parent == cellnode2.parent
            return cellnode1.parent
        #if both cells are roots there is no MRCA
        elseif isnothing(cellnode1.parent) && isnothing(cellnode2.parent)
            return nothing
        #if cellnode1 is parent of cellnode2 return cellnode1
        elseif cellnode1 == cellnode2.parent
            return cellnode1
        #if none of these are satisfied keep looking 
        else
            cellnode2 = cellnode2.parent
        end
    end
end

function findMRCA(cellnodes::Vector)
    cellnodes = copy(cellnodes)
    cellnode1 = pop!(cellnodes)
    while length(cellnodes) > 0
        cellnode2 = pop!(cellnodes)
        cellnode1 = findMRCA(cellnode1, cellnode2)
    end
    return cellnode1
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