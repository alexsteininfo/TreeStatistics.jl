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


## Things we need to define 
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
AbstractTrees.SiblingLinks(::Type{<:BinaryNode}) = AbstractTrees.StoredSiblings()
AbstractTrees.parent(n::BinaryNode) = n.parent
AbstractTrees.NodeType(::Type{<:BinaryNode{T}}) where {T} = HasNodeType()
AbstractTrees.nodetype(::Type{<:BinaryNode{T}}) where {T} = BinaryNode{T}

#ensure all nodes in tree are of the same type
Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()

AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

# function AbstractTrees.printnode(io::IO, node::BinaryNode)
#     print(io, "$(node.data.id), $(node.data.mutations) ($(popsize(node)))")
#     node.data.alive || print(io, " X")
# end

function popsize(root::BinaryNode)
    N = 0
    for l in Leaves(root)
        if alive(nodevalue(l))
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

AbstractTrees.getroot(nodevec::Vector{BinaryNode{T}}) where T = collect(Set(getroot.(nodevec)))

function getsingleroot(nodevec::Vector{BinaryNode{T}}) where T 
    roots = AbstractTrees.getroot(nodevec)
    if length(roots) == 1
        return roots[1]
    else
        return nothing
    end
end

function findMRCA(cellnode1, cellnode2)
    if cellnode1.data.id > cellnode2.data.id
        cellnode1, cellnode2 = cellnode2, cellnode1
    end
    if cellnode1 == cellnode2
        return cellnode1
    elseif cellnode1.parent == cellnode2.parent
        return cellnode1.parent
    elseif isdefined(cellnode2, :parent)
         return findMRCA(cellnode1, cellnode2.parent)
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