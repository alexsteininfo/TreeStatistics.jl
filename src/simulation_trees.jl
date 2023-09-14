popsize_exceeded(popsize, input::BranchingInput) = popsize > input.Nmax
popsize_exceeded(popsize, input::SinglelevelInput) = false

function simulate!(treemodule::TreeModule, input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, t0=nothing, tmax=nothing)
    branchingprocess!(
        treemodule,
        input.birthrate, 
        input.deathrate, 
        input.Nmax, 
        input.μ, 
        input.mutationdist, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        rng;
        timefunc,
        t0
    )
    return treemodule
end

function simulate!(treemodule::TreeModule, input::MoranInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, t0=nothing, tmax=nothing)
    moranprocess!(
        treemodule,
        input.moranrate, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        input.μ, 
        input.mutationdist, 
        rng;
        timefunc,
        t0,
        moranincludeself=input.moranincludeself
    )
    return treemodule
end

function simulate!(treemodule::TreeModule, input::BranchingMoranInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, t0=nothing, tmax=nothing)
    
    if length(treemodule) < input.Nmax
        branchingprocess!(
            treemodule,
            input.birthrate, 
            input.deathrate, 
            input.Nmax, 
            input.μ, 
            input.mutationdist, 
            isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
            rng; 
            timefunc,
            t0
        )
        t0 = age(treemodule)
    end

    moranprocess!(
        treemodule,
        input.moranrate, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        input.μ, 
        input.mutationdist, 
        rng;
        timefunc,
        t0,
        moranincludeself=input.moranincludeself
    )
    return treemodule
end

"""
    branchingprocess!(treemodule::TreeModule, birthrate, deathrate, Nmax, μ, mutationdist, 
        tmax, rng::AbstractRNG; timefunc=exptime)

Simulate a population of cells, defined by `treemodule` that grows by a branching process.

"""
function branchingprocess!(treemodule::TreeModule, birthrate, deathrate, Nmax, μ, mutationdist, 
    tmax, rng::AbstractRNG; timefunc=exptime, t0=nothing)

    # set initial time, population size and next cell ID
    t = !isnothing(t0) ? t0 : maximum(cellnode.data.birthtime for cellnode in treemodule.cells)

    N = length(treemodule.cells)
    nextID = maximum(cellnode.data.id for cellnode in treemodule.cells)

    while N < Nmax && N > 0
        Δt = timefunc(rng, N * (birthrate + deathrate))
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt
        _, N, nextID = 
            branchingupdate!(treemodule, birthrate, deathrate, N, t, nextID, μ, mutationdist, rng)
    end
    #add final mutations to all alive cells if mutations are time dependent
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        add_mutations!(treemodule, t, mutationdist, μ, rng)
    end

    return treemodule

end

"""
    branchingupdate!(treemodule::TreeModule, birthrate, deathrate, N, t, nextID, μ, mutationdist, rng; 
    timefunc=exptime)
    
Single update step of branching process.
"""
function branchingupdate!(treemodule::TreeModule, birthrate, deathrate, N, t, nextID, μ, 
    mutationdist, rng)

    #pick a random cell and randomly select its fate (birth or death) with probability 
    #proportional to birth and death rates
    randcellidx = rand(rng, 1:N) 
    r = rand(rng) 
    if r < birthrate / (birthrate + deathrate)
        _, nextID = 
            celldivision!(treemodule, randcellidx, t, nextID, μ, mutationdist, rng)
        N += 1
        updatetime!(treemodule, t)

    else
        celldeath!(treemodule, randcellidx, t, μ, mutationdist, rng)
        N -= 1
        updatetime!(treemodule, t)

    end
    return treemodule, N, nextID
end

"""
    moranprocess!(treemodule::TreeModule, moranrate, tmax, μ, mutationdist, rng; 
        N=length(treemodule), timefunc=exptime)

Simulate a population of cells in `treemodule` with Moran process dynamics.
"""
function moranprocess!(treemodule::TreeModule, moranrate, tmax, μ, mutationdist, rng; 
    N=length(treemodule), timefunc=exptime, t0=nothing, moranincludeself=true)

    # set initial time and next cell ID
    t = !isnothing(t0) ? t0 : maximum(cellnode.data.birthtime for cellnode in treemodule.cells)
    nextID = maximum(cellnode.data.id for cellnode in treemodule.cells)

    while true
        Δt = timefunc(rng, N * moranrate)
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt
        _, nextID = 
            moranupdate!(treemodule, t, nextID, μ, mutationdist, rng; N, timefunc, moranincludeself)
    end
    #add final mutations to all alive cells if mutations are time dependent
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        add_mutations!(treemodule, t, mutationdist, μ, rng)
    end

    return treemodule

end
"""
    moranupdate!(treemodule::TreeModule, t, nextID, μ, mutationdist, rng; 
        N=length(treemodule), timefunc=timefunc)

Single update step of Moran process.
"""
function moranupdate!(treemodule::TreeModule, t, nextID, μ, mutationdist, rng; 
    N=length(treemodule), timefunc=timefunc, moranincludeself=true)

    #pick a cell to divide and a cell to die
    dividecellidx = rand(rng, 1:N) 
    deadcellidx = 
        if moranincludeself 
            deadcellidx = rand(rng, 1:N)
            #if dead cell and divide cell are the same kill one of the offspring
            deadcellidx = deadcellidx == dividecellidx ? N + 1 : deadcellidx
        else
            #exclude dividecellidx
            deadcellidx = rand(rng, deleteat!(collect(1:N), dividecellidx))
        end

    _, nextID = celldivision!(treemodule, dividecellidx, t, nextID, μ, mutationdist, rng)
    celldeath!(treemodule, deadcellidx, t, μ, mutationdist, rng)
    updatetime!(treemodule, t)
    return treemodule, nextID
end

function asymmetricupdate!(treemodule::TreeModule, t, nextID, μ, mutationdist, rng; 
    N=length(treemodule), timefunc=timefunc)

    #pick a cell to divide
    dividecellidx = rand(rng, 1:N) 
    _, nextID = celldivision!(treemodule, dividecellidx, t, nextID, μ, mutationdist, rng; 
        nchildcells=1)
        updatetime!(cellmodule, t)
    return treemodule, nextID
end

"""
    add_mutations!(treemodule, t, mutationdist, μ, rng)

Add mutations to live cells in `treemodule`, that have occured between cell birth time and 
time `t`.

"""
function add_mutations!(treemodule, t, mutationdist, μ, rng)
    for cellnode in treemodule.cells
        Δt = t - cellnode.data.birthtime
        cellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    end
    return cells
end


"""
    celldivision!(alivecells::Vector{BinaryNode{T}}, parentcellidx, t, nextID, μ, 
    mutationdist, rng; nchildcells = 2) where T <: AbstractTreeCell

Remove parent cell from `alivecells` and add `nchildcells` new child cells (can be 1 or 2). 

If mutations are time-dependent, e.g. `mutationdist == poissontimedep`, add mutations to the
parent cell depending on the length of its lifetime. Otherwise assign mutations to each 
child cell.
"""
function celldivision!(treemodule::TreeModule{T}, parentcellidx, t, nextID, μ, 
    mutationdist, rng; nchildcells=2) where T <: AbstractTreeCell

    alivecells = treemodule.cells
    parentcellnode = alivecells[parentcellidx] #get parent cell node
    # deleteat!(alivecells, parentcellidx) #delete parent cell node from alivecells list

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
    alivecells[parentcellidx] = leftchild!(parentcellnode, childcell1)

    if nchildcells == 2
        childcell2 = T(
            id=nextID + 1,
            birthtime=t, 
            mutations=childcellmuts[2], 
            clonetype=parentcellnode.data.clonetype
        )
        #the other cell is added to the end of alivecells
        push!(alivecells, rightchild!(parentcellnode, childcell2)) 
    end
    return treemodule, nextID + nchildcells
end

"""
    celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellidx[, t, μ, mutationdist, rng]) 
        where T<: AbstractTreeCell

Remove dead cell from `alivecells` and remove references to it from tree. Add time dependent mutations (if applicable) to dying cell.
"""
function celldeath!(treemodule::TreeModule, deadcellidx, t=nothing, 
    μ=nothing, mutationdist=nothing, rng=nothing)

    alivecells = treemodule.cells
    #remove references to dead cell
    killcell!(alivecells, deadcellidx, t, μ, mutationdist, rng)
    #remove from alivecell vector
    deleteat!(alivecells, deadcellidx)
    return treemodule
end

# """
#     celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellidx::Vector[, t, μ, mutationdist, rng]) 
#         where T<: AbstractTreeCell
# """
# function celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellidx::Vector{Int64}, t=nothing, μ=nothing, 
#     mutationdist=nothing, rng=nothing) where T<:AbstractTreeCell
#     for id in deadcellidx
#         celldeath!(alivecells, id, t, μ, mutationdist, rng)
#     end
#     return alivecells
# end

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
    cellremoval!(cellmodule::TreeModule, deadcells)

Remove cells from module without killing them.
"""
function cellremoval!(cellmodule::TreeModule, deadcells)
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