mutable struct BinaryNode{T}
    data::T
    parent::BinaryNode{T}
    left::BinaryNode{T}
    right::BinaryNode{T}

    # Root constructor
    BinaryNode{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNode{T}(data, parent::BinaryNode{T}) where T = new{T}(data, parent)
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

function leftchild(data, parent::BinaryNode)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end
function rightchild(data, parent::BinaryNode)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end


## Enhancement of the "native" binary tree
# You might have the methods below even if you weren't trying to support AbstractTrees.

# Implement iteration over the immediate children of a node
function Base.iterate(node::BinaryNode)
    isdefined(node, :left) && return (node.left, false)
    isdefined(node, :right) && return (node.right, true)
    return nothing
end
function Base.iterate(node::BinaryNode, state::Bool)
    state && return nothing
    isdefined(node, :right) && return (node.right, true)
    return nothing
end
Base.IteratorSize(::Type{BinaryNode{T}}) where T = Base.SizeUnknown()
Base.eltype(::Type{BinaryNode{T}}) where T = BinaryNode{T}

## Things we need to define to leverage the native iterator over children
## for the purposes of AbstractTrees.
# Set the traits of this kind of tree
Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()
AbstractTrees.parentlinks(::Type{BinaryNode{T}}) where T = AbstractTrees.StoredParents()
AbstractTrees.siblinglinks(::Type{BinaryNode{T}}) where T = AbstractTrees.StoredSiblings()


function AbstractTrees.children(node::BinaryNode)
    if isdefined(node, :left)
        if isdefined(node, :right)
            return (node.left, node.right)
        end
        return (node.left,)
    end
    isdefined(node, :right) && return (node.right,)
    return ()
end

AbstractTrees.parent(node::BinaryNode) = isdefined(node, :parent) ? node.parent : nothing

function AbstractTrees.nextsibling(child::BinaryNode)
    isdefined(child, :parent) || return nothing
    p = child.parent
    if isdefined(p, :right)
        child === p.right && return nothing
        return p.right
    end
    return nothing
end



# We also need `pairs` to return something sensible.
# If you don't like integer keys, you could do, e.g.,
#   Base.pairs(node::BinaryNode) = BinaryNodePairs(node)
# and have its iteration return, e.g., `:left=>node.left` and `:right=>node.right` when defined.
# But the following is easy:
Base.pairs(node::BinaryNode) = enumerate(node)


AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

function AbstractTrees.printnode(io::IO, node::BinaryNode{SimpleCell})
    print("$(node.data.id), $(node.data.mutations) ($(popsize(node)))")
    node.data.alive || print(" X")
end

function popsize(root::BinaryNode)
    N = 0
    for l in Leaves(root)
        N += 1
    end
    return N
end

function popsize(root::BinaryNode, n)
    if n == 1
        return size(root)
    else
        return popsize(root.left, n-1), popsize(root.right, n-1)
    end
end

root(node::BinaryNode) = AbstractTrees.isroot(node) ? node : root(AbstractTrees.parent(node))
