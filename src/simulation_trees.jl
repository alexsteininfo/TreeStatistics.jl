popsize_exceeded(popsize, input::BranchingInput) = popsize > input.Nmax
popsize_exceeded(popsize, input::SinglelevelInput) = false

"""
    celldivision!(treemodule::TreeModule{T, S}, subclones, parentcellidx, t, nextID, μ, 
    mutationdist, rng; nchildcells=2) where {T <: AbstractTreeCell, S}

Remove parent cell from `alivecells` and add `nchildcells` new child cells (can be 1 or 2). 

If mutations are time-dependent, e.g. `mutationdist == poissontimedep`, add mutations to the
parent cell depending on the length of its lifetime. Otherwise assign mutations to each 
child cell.
"""
function celldivision!(treemodule::TreeModule{T, S}, subclones, parentcellid, t, nextID, μ, 
    mutationdist, rng; nchildcells=2) where {T <: AbstractTreeCell, S}

    alivecells = treemodule.cells
    parentcellnode = alivecells[parentcellid] #get parent cell node

    #if mutations are time dependent assign mutations to parent cell and give new cells
    #no initial mutations
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        childcellmuts = zeros(Int64, nchildcells)
        Δt = t - parentcellnode.data.birthtime
        parentcellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    #if mutations are not time dependent assign new child cell mutations
    else
        childcellmuts = [numbernewmutations(rng, mutationdist, μ) for _ in 1:nchildcells]
    end
    #create new child cells and add them to alivecells list
    childcell1 = T(
            id=nextID,
            birthtime=t, 
            mutations=childcellmuts[1], 
            clonetype=parentcellnode.data.clonetype
    )
    #one cell replaces the parent in alivecells
    alivecells[parentcellid] = leftchild!(parentcellnode, childcell1)

    if nchildcells == 2
        childcell2 = T(
            id=nextID + 1,
            birthtime=t, 
            mutations=childcellmuts[2], 
            clonetype=parentcellnode.data.clonetype
        )
        #the other cell is added to the end of alivecells
        push!(alivecells, rightchild!(parentcellnode, childcell2)) 
        subclones[parentcellnode.data.clonetype].size += 1
    end
    #adjust subclone sizes
    return treemodule, subclones, nextID + nchildcells
end

"""
    celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellidx[, t, μ, mutationdist, rng]) 
        where T<: AbstractTreeCell

Remove dead cell from `alivecells` and remove references to it from tree. Add time dependent mutations (if applicable) to dying cell.
"""
function celldeath!(treemodule::TreeModule, subclones::Vector{Subclone}, deadcellidx, t, 
    μ=nothing, mutationdist=nothing, rng=nothing)

    alivecells = treemodule.cells
    deadcellclonetype = alivecells[deadcellidx].data.clonetype
    #remove references to dead cell
    killcell!(alivecells, deadcellidx, t, μ, mutationdist, rng)
    #remove from alivecell vector
    deleteat!(alivecells, deadcellidx)
    #adjust subclone sizes
    subclones[deadcellclonetype].size -= 1
    return treemodule, subclones
end

"""
    killcell!(alivecells::Vector{BinaryNode{TreeCell}}, deadcellidx::Int64, t, μ, mutationdist, rng)

Kill `TreeCell` at index `deadcellidx` in `alivecells` by adding a left child node with `alive=false`. If 
`mutationdist ∈ [:poissontimedep, :fixedtimedep]` then add mutations to the dying cell.

"""
function killcell!(alivecells::TreeCellVector, deadcellidx::Int64, t, μ, mutationdist, rng)
    deadcellnode = alivecells[deadcellidx]
    #if mutations are time dependent, add the number accumulated by the cell
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep
        Δt = t - deadcellnode.data.birthtime
        deadcellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    end
    #add leftchild node containing a dead cell
    leftchild!(deadcellnode, TreeCell(deadcellnode.data.id, false, t, 0, deadcellnode.data.clonetype))
end

"""
    killcell!(alivecells::Vector{BinaryNode{SimpleTreeCell}}, deadcellidx::Int64)

Kill `SimpleTreeCell` at index `deadcellidx` in `alivecells` by removing all references to it from the tree.
"""
function killcell!(alivecells::SimpleTreeCellVector, deadcellidx::Int64, args...)
    prune_tree!(alivecells[deadcellidx])
    return alivecells
end

"""
    killcells!(alivecells, deadcellvector, args...)

Kill 
"""
function killcells!(alivecells, deadcellvector, args...)
    for deadcellidx in deadcellvector
        killcell!(alivecells, deadcellidx, args...)
    end
    return alivecells
end

function killallcells!(alivecells, args...)
    killcells!(alivecells, 1:lastindex(alivecells), args...)
    return alivecells
end

"""
    cellremoval!(module, deadcells)

Remove cells from module without killing them.
"""
function cellremoval!(cellmodule, deadcells)
    deleteat!(cellmodule.cells, sort!(deadcells))
    return cellmodule.cells
end

    

"""
    prune_tree!(cellnode)

Remove `cellnode` from tree, and remove any node that has no children after its 
removal.

"""
function prune_tree!(cellnode)
    while true
        parent = cellnode.parent
        if isnothing(parent)
            return
        else
            cellnode.parent = nothing
            if parent.left == cellnode
                parent.left = nothing
            elseif parent.right == cellnode
                parent.right = nothing
            else
                error("dead cell is neither left nor right child of parent")
            end
            #if parent cell now has no children, it becomes cellnode to be removed
            if isnothing(parent.left) && isnothing(parent.right)
                cellnode = parent
            else
                return
            end
        end
    end
end


"""
    changemutations!(root::BinaryNode, μ, mutationdist, tmax, rng, clonalmutations=0)

Take the phylogeny, defined by `root` and assign new mutations according the the
`mutationdist` and other parameters.

"""
function changemutations!(root::BinaryNode, μ, mutationdist, tmax, rng, clonalmutations=0)
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        for cellnode in PreOrderDFS(root)
            Δt = celllifetime(cellnode, tmax)
            cellnode.data.mutations = numbernewmutations(rng, mutationdist, μ, Δt=Δt)
        end
    else
        for cellnode in PreOrderDFS(root)
            cellnode.data.mutations = numbernewmutations(rng, mutationdist, μ)
        end
    end
    root.data.mutations += clonalmutations
end

"""
    endtime(cellnode::BinaryNode)

Return the time at which the cell divided (or died if it is a `TreeCell`.) If the cell is 
still alive, return `nothing`.
"""
function endtime(cellnode::BinaryNode)
    if haschildren(cellnode)
        return cellnode.left.data.birthtime
    else
        return nothing
    end
end 

"""
    celllifetime(cellnode::BinaryNode, [tmax])

Comuptes the lifetime of a given cell. If it hasn't yet died or divided, end of lifetime is
given by tmax (defaults to the age of the population).

"""
function celllifetime(cellnode::BinaryNode, tmax=nothing)
    if haschildren(cellnode)
        return cellnode.left.data.birthtime - cellnode.data.birthtime
    else
        tmax = isnothing(tmax) ? age(getroot(cellnode)) : tmax
        return tmax - cellnode.data.birthtime

    end
end

"""
    celllifetimes(root; excludeliving=true)

Computes the lifetime of each cell in the phylogeny, excluding currently alive cells by
default.    
"""
function celllifetimes(root; excludeliving=true)
    lifetimes = Float64[]
    if excludeliving
        for cellnode in PreOrderDFS(root)
            if haschildren(cellnode)
                push!(lifetimes, cellnode.left.data.birthtime - cellnode.data.birthtime)
            end
        end
    else
        popage = age(root)
        for cellnode in PreOrderDFS(root)
            push!(lifetimes, celllifetime(cellnode, popage))
        end
    end
    return lifetimes
end

"""
    age(root::BinaryNode)

Compute the age of the population, given by the time of the most recent cell division.
"""
function age(root::BinaryNode)
    age = 0
    for cellnode in Leaves(root)
        if cellnode.data.birthtime > age
            age = cellnode.data.birthtime
        end
    end
    return age
end

getalivecells(root::BinaryNode) = 
    [cellnode for cellnode in Leaves(root) if isalive(cellnode.data)]

getalivecells(roots::Vector{BinaryNode{T}}) where T = 
    [cellnode for root in roots for cellnode in Leaves(root) if isalive(cellnode.data)]

popsize(root::BinaryNode{SimpleTreeCell}) = treebreadth(root)
popsize(roots::Vector) = sum(popsize(root) for root in roots)

isalive(cellnode::BinaryNode{T}) where T = isalive(cellnode.data)
isalive(cell::TreeCell) = cell.alive
isalive(cell::SimpleTreeCell) = true
isalive(::Nothing) = false

id(cellnode::BinaryNode{<:AbstractTreeCell}) = cellnode.data.id


"""
    asroot!(node)

Transform `node` into a root by setting `parent` field to nothing. Return `node` and the 
original `parent` node.
"""
function asroot!(node)
    parent = node.parent
    node.parent = nothing
    return node, parent
end

"""
    cell_subset_size(node, cells)

Calculate the number of cell nodes (leaves of the tree of which `node` is the root) that are
both alive and present in the list `cells`.
"""
function cell_subset_size(node, cells)
    if node in cells
        if isalive(node) return 1 else 0 end
    end
    #To properly iterate over leaves we need to make node a true "root" (i.e. set its parent
    #field to nothing)
    root, parent = asroot!(node)
    count = 0
    for leaf in Leaves(root)
        if leaf in cells && isalive(leaf)
            count += 1
        end
    end
    root.parent = parent #reset node so that it is unchanged
    return count
end

getnextID(cells::TreeCellVector) = maximum(cellnode.data.id for cellnode in cells)