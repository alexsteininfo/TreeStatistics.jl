"""
    run1simulation_tree(input::BranchingInput, [rng::AbstractRNG]; timefunc=exptime)

    Simulate a population of cells that grows by a branching process.

    Take simulation parameters from `input` and return a BinaryNode that is the root of a 
    phylogeny tree. By default if population goes extinct, restart simulation. 

"""
function run1simulation_tree(input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, returnextinct=false)
    
    while true
        alivecells, root = initialize_tree(input.clonalmutations)
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
            return alivecells, root
        end
    end
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
    branchingupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, N, t, nextID, μ, mutationdist, 
        rng; timefunc=exptime)
    
    Single update step of branching process.
"""
function branchingupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, N, t, 
    nextID, μ, mutationdist, rng; timefunc=exptime)

    randcellidx = rand(rng, 1:N)
    r = rand(rng)
    if r < b/(b+d)
        N, nextID = 
            celldivision!(alivecells, randcellidx, N, t, nextID, μ, mutationdist, rng)
    else
        N = celldeath!(alivecells, randcellidx, N, t, μ, mutationdist, rng)
    end
    return N, nextID
end

"""
    celldivision!(alivecells::Vector{BinaryNode{SimpleCell}}, parentcellidx, N, t, nextID, μ, 
        mutationdist, rng)
Remove parent cell from `alivecells` vector and add two new child cells. 

If mutations are time-dependent, e.g. `mutationdist == poissontimedep`, add mutations to the
parent cell depending on the length of its lifetime. Otherwise assign mutations to each 
child cell.
"""

function celldivision!(alivecells::Vector{BinaryNode{SimpleCell}}, parentcellidx, N, t, nextID, μ, 
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

    return N + 1, nextID + 2

end


"""
    celldeath!(alivecells::Vector{BinaryNode{SimpleCell}}, deadcellidx, N, t, [μ, 
        mutationdist, rng])

Remove dead cell from `alivecells` vector and add a new cell with `alive=false` as the left
child of the dead cell. Add time dependent mutations (if applicable) to dying cell.
"""
function celldeath!(alivecells::Vector{BinaryNode{SimpleCell}}, deadcellidx, N, t, μ=nothing, 
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
    return N - 1   
end

"""
    initialize_tree(clonalmutations, [N])
Initialize tree with `N` cells (defaults to 1) and return vector of alive cells.
"""
function initialize_tree(clonalmutations, N=1)
    if N == 1
        initialcell = SimpleCell(1, true, 0.0, clonalmutations, 1)
        root = BinaryNode{SimpleCell}(initialcell)
        alivecells = [root]
        return alivecells, root
    else 
        initialcells = [SimpleCell(1, true, 0.0, clonalmutations, 1) for i in 1:N]
        roots = [BinaryNode{SimpleCell}(cell) for cell in initialcells]
        alivecells = copy(roots)
        return alivecells, roots
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