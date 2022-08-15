
"""
    run1simulation_tree(input::BranchingInput, [rng::AbstractRNG]; timefunc=exptime)

    Simulate a population of cells that grows by a branching process.

    Take simulation parameters from `input` and return `alivecells` (a vector of type 
    BinaryNode corresponding to the alive cells at the end of the simulation) and `root` a 
    single BinaryNode that is the root of a phylogeny tree. By default if population goes 
    extinct, restart simulation. 

"""
function run1simulation_tree(::Type{T}, input::SimulationInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, returnextinct=false) where T <: AbstractTreeCell
    while true
        alivecells = initialize(T, input, rng)
        alivecells = run1simulation_tree(alivecells, input, rng; timefunc)
        (length(alivecells) > 0 || returnextinct) && return alivecells
    end

end

function run1simulation_tree(cell::AbstractTreeCell, input::SimulationInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime)

    return run1simulation_tree([BinaryNode(cell)], input, rng; timefunc)

end

function run1simulation_tree(alivecells::Vector{BinaryNode{T}}, input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime) where T <: AbstractTreeCell
    
    alivecells =
        branchingprocess!(
            alivecells,
            input.b, 
            input.d, 
            input.Nmax, 
            input.μ, 
            input.mutationdist, 
            input.tmax, 
            rng, 
            timefunc=timefunc
        )
    return alivecells
end

"""
    run1simulation_tree(input::MoranInput, [rng::AbstractRNG]; timefunc=exptime)

    Simulate a population of cells that evolves by a Moran process.

    Take simulation parameters from `input` and return `alivecells` (a vector of type 
    BinaryNode corresponding to the alive cells at the end of the simulation) and `roots` (a 
    vector of BinaryNodes that are the roots of the phylogeny).
"""
function run1simulation_tree(alivecells::Vector{BinaryNode{T}}, input::MoranInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime) where T <: AbstractTreeCell
    
    alivecells =
        moranprocess!(
            alivecells,
            input.bdrate, 
            input.tmax, 
            input.μ, 
            input.mutationdist, 
            rng, 
            timefunc=timefunc
        )
    return alivecells
end

"""
    branchingprocess!(alivecells::Vector{BinaryNode{TreeCell}}, b, d, Nmax, μ, mutationdist, tmax,
        rng; timefunc=exptime)

    Simulate a population of cells that grows by a branching process, starting with the 
    population of cells given in the alivecells vector.

"""
function branchingprocess!(alivecells::Vector{BinaryNode{T}}, b, d, Nmax, μ, mutationdist, 
    tmax, rng::AbstractRNG; timefunc=exptime) where T <: AbstractTreeCell

    # set initial time, population size and next cell ID
    t = maximum(cellnode.data.birthtime for cellnode in alivecells)
    N = length(alivecells)
    nextID = maximum(cellnode.data.id for cellnode in alivecells)

    while N < Nmax && N > 0
        Δt = timefunc(rng, N * (b + d))
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt
        _, N, nextID = 
            branchingupdate!(alivecells, b, d, N, t, nextID, μ, mutationdist, rng, 
                timefunc=timefunc)
    end
    #add final mutations to all alive cells if mutations are time dependent
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        add_mutations!(alivecells, t, mutationdist, μ, rng)
    end

    return alivecells

end

"""
    branchingupdate!(alivecells::Vector{BinaryNode{TreeCell}}, b, d, N, t, nextID, μ, mutationdist, 
        rng; timefunc=exptime)
    
    Single update step of branching process.
"""
function branchingupdate!(alivecells::Vector{BinaryNode{T}}, b, d, N, t, 
    nextID, μ, mutationdist, rng; timefunc=exptime) where T <: AbstractTreeCell

    #pick a random cell and randomly select its fate (birth or death) with probability 
    #proportional to birth and death rates
    randcellidx = rand(rng, 1:N) 
    r = rand(rng) 
    if r < b/(b+d)
        _, nextID = 
            celldivision!(alivecells, randcellidx, t, nextID, μ, mutationdist, rng)
        N += 1
    else
        celldeath!(alivecells, randcellidx, t, μ, mutationdist, rng)
        N -= 1
    end
    return alivecells, N, nextID
end

"""
    moranprocess!(alivecells::Vector{BinaryNode{TreeCell}}, bdrate, tmax, μ, mutationdist, 
        rng::AbstractRNG; timefunc=exptime)

"""
function moranprocess!(alivecells::Vector{BinaryNode{T}}, bdrate, tmax, μ, mutationdist, 
    rng::AbstractRNG; N=length(alivecells), timefunc=exptime) where T <: AbstractTreeCell

    # set initial time and next cell ID
    t = maximum(cellnode.data.birthtime for cellnode in alivecells)
    nextID = maximum(cellnode.data.id for cellnode in alivecells)

    while true
        Δt = timefunc(rng, N * bdrate)
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt
        _, nextID = 
            moranupdate!(alivecells, t, nextID, μ, mutationdist, rng; N, timefunc)
    end
    #add final mutations to all alive cells if mutations are time dependent
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        add_mutations!(alivecells, t, mutationdist, μ, rng)
    end

    return alivecells

end
"""
    moranupdate!(alivecells::Vector{BinaryNode{TreeCell}}, t, nextID, μ, 
        mutationdist, rng, timefunc=timefunc)
"""
function moranupdate!(alivecells::Vector{BinaryNode{T}}, t, nextID, μ, 
    mutationdist, rng; N=length(alivecells), timefunc=timefunc) where T <: AbstractTreeCell

    #pick a cell to divide and a cell to die
    dividecellidx = rand(rng, 1:N) 
    deadcellidx = rand(rng, 1:N) 

    #if dead cell and divide cell are the same kill one of the new offspring after division
    if deadcellidx == dividecellidx 
        deadcellidx = N
    end

    _, nextID = celldivision!(alivecells, dividecellidx, t, nextID, μ, mutationdist, rng)
    celldeath!(alivecells, deadcellidx, t, μ, mutationdist, rng)
    return alivecells, nextID
end

"""
    add_mutations(cells, t)

Add mutations to vector of `cells`, that have occured between cell birth time and time t.

"""
function add_mutations!(cells, t, mutationdist, μ, rng)
    for cellnode in cells
        Δt = t - cellnode.data.birthtime
        cellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    end
    return cells
end


"""
    celldivision!(alivecells::Vector{BinaryNode{TreeCell}}, parentcellidx, N, t, nextID, μ, 
        mutationdist, rng)
Remove parent cell from `alivecells` vector and add two new child cells. 

If mutations are time-dependent, e.g. `mutationdist == poissontimedep`, add mutations to the
parent cell depending on the length of its lifetime. Otherwise assign mutations to each 
child cell.
"""

function celldivision!(alivecells::Vector{BinaryNode{T}}, parentcellidx, t, nextID, μ, 
    mutationdist, rng) where T <: AbstractTreeCell

    parentcellnode = alivecells[parentcellidx] #get parent cell node
    deleteat!(alivecells, parentcellidx) #delete parent cell node from alivecells list

    #if mutations are time dependent assign mutations to parent cell and give new cells
    #no initial mutations
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        childcellmuts1, childcellmuts2 = 0, 0
        Δt = t - parentcellnode.data.birthtime
        parentcellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    #if mutations are not time dependent assign new child cell mutations
    else
        childcellmuts1 = numbernewmutations(rng, mutationdist, μ)
        childcellmuts2 = numbernewmutations(rng, mutationdist, μ)
    end
    #create new child cells and add them to alivecells list
    childcell1 = T(id=nextID, birthtime=t, mutations=childcellmuts1, 
        clonetype=parentcellnode.data.clonetype)    
    childcell2 = T(id=nextID + 1, birthtime=t, mutations=childcellmuts2, 
        clonetype=parentcellnode.data.clonetype)
    push!(alivecells, leftchild!(parentcellnode, childcell1))
    push!(alivecells, rightchild!(parentcellnode, childcell2)) 

    return alivecells, nextID + 2

end

function celldivision!(moduletracker::TreeModule{T}, parentcellidx, t, nextID, μ, 
    mutationdist, rng) where T <: AbstractTreeCell

    return celldivision!(moduletracker.cells, parentcellidx, t, nextID, μ, mutationdist, rng)

end
"""
    celldeath!(alivecells::Vector{BinaryNode{TreeCell}}, deadcellidx, t, μ=nothing, 
        mutationdist=nothing, rng=nothing)

Remove dead cell from `alivecells` vector and add a new cell with `alive=false` as the left
child of the dead cell. Add time dependent mutations (if applicable) to dying cell.
"""
function celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellidx::Int64, t, μ=nothing, 
    mutationdist=nothing, rng=nothing) where T

    killcell!(alivecells, deadcellidx, t, μ, mutationdist, rng)
    #remove from alivecell vector
    deleteat!(alivecells, deadcellidx)

    return alivecells
end

function celldeath!(alivecells::Vector{BinaryNode{SimpleTreeCell}}, deadcellidx::Int64)

    #remove references to dead cell
    killcell!(alivecells::Vector{BinaryNode{SimpleTreeCell}}, deadcellidx::Int64)
    #remove from alivecell vector
    deleteat!(alivecells, deadcellidx)
    return alivecells
end

"""
    celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellvector::Vector{Int64}, args...) where T<:AbstractTreeCell

TBW
"""
function celldeath!(alivecells::Vector{BinaryNode{T}}, deadcellvector::Vector{Int64}, args...) where T<:AbstractTreeCell
    for id in deadcellvector
        celldeath!(alivecells, id, args...)
    end
    return alivecells
end

"""
    celldeath!(moduletracker::TreeModule, args...)

"""
function celldeath!(moduletracker::TreeModule, args...)
    celldeath!(moduletracker.cells, args...)
    return moduletracker
end

function killcell!(alivecells::Vector{BinaryNode{TreeCell}}, deadcellidx::Int64, t, μ, mutationdist, rng)
    deadcellnode = alivecells[deadcellidx]
    #if mutations are time dependent, add the number accumulated by the cell
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep
        Δt = t - deadcellnode.data.birthtime
        deadcellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    end
    #add leftchild node containing a dead cell
    leftchild!(deadcellnode, TreeCell(deadcellnode.data.id, false, t, 0, deadcellnode.data.clonetype))
end

function killcell!(alivecells::Vector{BinaryNode{SimpleTreeCell}}, deadcellidx::Int64, args...)
    prune_tree!(alivecells[deadcellidx])
    return alivecells
end

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
    cellremoval!(moduletracker::TreeModule, deadcells)

Remove cells from module without killing them.
"""
function cellremoval!(moduletracker::TreeModule, deadcells)
    deleteat!(moduletracker.cells, deadcells)
    return moduletracker.cells
end

    

"""
    prune_tree!(cellnode)

Remove `cellnode` from tree, and recursively remove any node that has no children after its 
removal.

"""
function prune_tree!(cellnode)
    if !isnothing(cellnode.parent)
        parent = cellnode.parent
        cellnode.parent = nothing
        if parent.left == cellnode
            parent.left = nothing
        elseif parent.right == cellnode
            parent.right = nothing
        else
            error("dead cell is neither left nor right child of parent")
        end
        if isnothing(parent.left) && isnothing(parent.right)
            prune_tree!(parent)
        end
    end
end

"""
    initialize(clonalmutations, [N])
Initialize tree with `N` cells (defaults to 1) and return vector of alive cells.
"""
function initialize(::Type{T}, initialmutations=0, N=1) where T <:AbstractTreeCell
    alivecells = map(
        id -> BinaryNode{T}(T(id=id, mutations=initialmutations)), 
        1:N
    )
    return alivecells
end

function initialize(::Type{T}, input::MoranInput, rng) where T <: AbstractTreeCell
    if input.mutationdist == :fixedtimedep || input.mutationdist == :poissontimedep    
        alivecells = map(
            id -> BinaryNode{T}(T(id=id, mutations=0)), 
            1:input.N
        )
    else 
        alivecells = map(
            id -> BinaryNode{T}(T(
                id=id, 
                mutations=numbernewmutations(rng, input.mutationdist, input.μ)
            )), 
            1:input.N
        )
    end
    return alivecells
end

function initialize(::Type{T}, input::Union{BranchingInput, BranchingMoranInput, MultilevelInput}, rng) where T <: AbstractTreeCell
    if input.mutationdist == :fixedtimedep || input.mutationdist == :poissontimedep    
        alivecells = [BinaryNode{T}(T(id=1, mutations=0))]
    else 
        alivecells = [BinaryNode{T}(T(
            id=1, 
            mutations=numbernewmutations(rng, input.mutationdist, input.μ)
        ))]          
    end
    return alivecells
end

function initializesim_from_cells(cells::Vector{BinaryNode{T}}, subclones::Vector{CloneTracker}, id, parentid; inittime=0.0) where T<: AbstractTreeCell
    return initializesim_from_cells(
        TreeModule{T}, 
        cells, 
        subclones,
        id, 
        parentid; 
        inittime
    )
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

Return the time at which the cell divided or died. If the cell is still alive, return 
`nothing`.
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
    [cellnode for cellnode in Leaves(root) if alive(cellnode.data)]

getalivecells(roots::Vector{BinaryNode{T}}) where T = 
    [cellnode for root in roots for cellnode in Leaves(root) if alive(cellnode.data)]

popsize(root::BinaryNode{SimpleTreeCell}) = treebreadth(root)
popsize(roots::Vector) = sum(popsize(root) for root in roots)

alive(cell::TreeCell) = cell.alive
alive(cell::SimpleTreeCell) = true

id(cellnode::BinaryNode{<:AbstractTreeCell}) = cellnode.data.id