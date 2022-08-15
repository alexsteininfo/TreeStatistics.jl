"""
    run1simulation_tree(input::BranchingInput, [rng::AbstractRNG]; timefunc=exptime)

    Simulate a population of cells that grows by a branching process.

    Take simulation parameters from `input` and return `alivecells` (a vector of type 
    BinaryNode corresponding to the alive cells at the end of the simulation) and `root` a 
    single BinaryNode that is the root of a phylogeny tree. By default if population goes 
    extinct, restart simulation. 

"""
function run1simulation_tree(input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, returnextinct=false)
    
    while true
        alivecells = initialize_tree(input.clonalmutations)
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
        if length(alivecells) > 0 || returnextinct
            return alivecells, getroot(alivecells[1])
        end
    end
end

"""
    run1simulation_tree(input::MoranInput, [rng::AbstractRNG]; timefunc=exptime)

    Simulate a population of cells that evolves by a Moran process.

    Take simulation parameters from `input` and return `alivecells` (a vector of type 
    BinaryNode corresponding to the alive cells at the end of the simulation) and `root` (a 
    vector of BinaryNodes that are the roots of the phylogeny).
"""
function run1simulation_tree(input::MoranInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime)
    
    alivecells = initialize_tree(input.clonalmutations)
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
end

"""
    branchingprocess!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, Nmax, μ, mutationdist, tmax,
        rng; timefunc=exptime)

    Simulate a population of cells that grows by a branching process, starting with the 
    population of cells given in the alivecells vector.

"""
function branchingprocess!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, Nmax, μ, mutationdist, 
    tmax, rng::AbstractRNG; timefunc=exptime)

    # set initial time, population size and next cell ID
    t = maximum(cellnode.data.birthtime for cellnode in alivecells)
    N = length(alivecells)
    nextID = maximum(cellnode.data.id for cellnode in alivecells)

    while N < Nmax && N > 0
        Δt = timefunc(rng, N * (b + d))
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt
        N, nextID = 
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
    branchingupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, N, t, nextID, μ, mutationdist, 
        rng; timefunc=exptime)
    
    Single update step of branching process.
"""
function branchingupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, N, t, 
    nextID, μ, mutationdist, rng; timefunc=exptime)

    #pick a random cell and randomly select its fate (birth or death) with probability 
    #proportional to birth and death rates
    randcellidx = rand(rng, 1:N) 
    r = rand(rng) 
    if r < b/(b+d)
        nextID = 
            celldivision!(alivecells, randcellidx, t, nextID, μ, mutationdist, rng)
        N += 1
    else
        celldeath!(alivecells, randcellidx, t, μ, mutationdist, rng)
        N -= 1
    end
    return N, nextID
end

"""
    moranprocess!(alivecells::Vector{BinaryNode{SimpleCell}}, bdrate, tmax, μ, mutationdist, 
        rng::AbstractRNG; timefunc=exptime)

"""
function moranprocess!(alivecells::Vector{BinaryNode{SimpleCell}}, bdrate, tmax, μ, mutationdist, 
    rng::AbstractRNG; N=length(alivecells), timefunc=exptime)

    # set initial time and next cell ID
    t = maximum(cellnode.data.birthtime for cellnode in alivecells)
    nextID = maximum(cellnode.data.id for cellnode in alivecells)

    while true
        Δt = timefunc(rng, N * bdrate)
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt
        nextID = 
            moranupdate!(alivecells, t, nextID, μ, mutationdist, rng; N, timefunc)
    end
    #add final mutations to all alive cells if mutations are time dependent
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
        add_mutations!(alivecells, t, mutationdist, μ, rng)
    end

    return alivecells

end
"""
    moranupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, t, nextID, μ, 
        mutationdist, rng, timefunc=timefunc)
"""
function moranupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, t, nextID, μ, 
    mutationdist, rng; N=length(alivecells), timefunc=timefunc)

    #pick a cell to divide and a cell to die
    dividecellidx = rand(rng, 1:N) 
    deadcellidx = rand(rng, 1:N) 

    #if dead cell and divide cell are the same kill one of the new offspring after division
    if deadcellidx == dividecellidx 
        deadcellidx = N
    end

    nextID = celldivision!(alivecells, dividecellidx, t, nextID, μ, mutationdist, rng)
    celldeath!(alivecells, deadcellidx, t, μ, mutationdist, rng)
    return nextID
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
end


"""
    celldivision!(alivecells::Vector{BinaryNode{SimpleCell}}, parentcellidx, N, t, nextID, μ, 
        mutationdist, rng)
Remove parent cell from `alivecells` vector and add two new child cells. 

If mutations are time-dependent, e.g. `mutationdist == poissontimedep`, add mutations to the
parent cell depending on the length of its lifetime. Otherwise assign mutations to each 
child cell.
"""

function celldivision!(alivecells::Vector{BinaryNode{SimpleCell}}, parentcellidx, t, nextID, μ, 
    mutationdist, rng)

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
    childcell1 = SimpleCell(nextID, true, t, childcellmuts1, parentcellnode.data.clonetype)    
    childcell2 = SimpleCell(nextID + 1, true, t, childcellmuts2, parentcellnode.data.clonetype)
    push!(alivecells, leftchild(childcell1, parentcellnode))
    push!(alivecells, rightchild(childcell2, parentcellnode)) 

    return nextID + 2

end


"""
    celldeath!(alivecells::Vector{BinaryNode{SimpleCell}}, deadcellidx, N, t, [μ, 
        mutationdist, rng])

Remove dead cell from `alivecells` vector and add a new cell with `alive=false` as the left
child of the dead cell. Add time dependent mutations (if applicable) to dying cell.
"""
function celldeath!(alivecells::Vector{BinaryNode{SimpleCell}}, deadcellidx, t, μ=nothing, 
    mutationdist=nothing, rng=nothing)

    deadcellnode = alivecells[deadcellidx]
    #if mutations are time dependent, add the number accumulated by the cell
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep
        Δt = t - deadcellnode.data.birthtime
        deadcellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    end

    #remove from alivecell vector and add leftchild node containing a dead cell
    deleteat!(alivecells, deadcellidx)
    leftchild(SimpleCell(deadcellnode.data.id, false, t, 0, deadcellnode.data.clonetype), deadcellnode)
    return
end

"""
    initialize_tree(clonalmutations, [N])
Initialize tree with `N` cells (defaults to 1) and return vector of alive cells.
"""
function initialize_tree(clonalmutations, N=1)
    alivecells = map(
        id -> BinaryNode{SimpleCell}(SimpleCell(id=id, mutations=clonalmutations)), 
        1:N
    )
    return alivecells
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
    [cellnode for cellnode in Leaves(root) if cellnode.data.alive]