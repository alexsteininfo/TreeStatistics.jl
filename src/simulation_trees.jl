function run1simulation_tree(input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime)
    
    alivecells, root = 
        generate_branching_tree(
            input.b, 
            input.d, 
            input.Nmax, 
            input.μ, 
            input.mutationdist, 
            input.clonalmutations, 
            input.tmax, 
            rng, 
            timefunc=timefunc
        )
    return alivecells, root
end

function generate_branching_tree(b, d, Nmax, μ, mutationdist, clonalmutations, tmax, rng::AbstractRNG; timefunc=exptime)
    alivecells, root, N, nextID = initialize_branching_tree(clonalmutations, 1)
    t = 0.0
    while N < Nmax && N > 0
        t += timefunc(rng, N * (b + d))
        if t > tmax
            break
        end
        N, nextID = 
            branchingupdate!(alivecells, b, d, N, t, nextID, μ, mutationdist, rng, timefunc=timefunc)
    end
    #add final mutations to all Cells
    for cellnode in alivecells
        Δt = t - cellnode.data.birthtime
        cellnode.data.mutations += numbernewmutations(rng, mutationdist, μ, Δt=Δt)
    end
    return alivecells, root

end

#currently for neutral evolution only, add selection later
function branchingupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, N, t, 
    nextID, μ, mutationdist, rng; timefunc=exptime)

    randcellidx = rand(rng, 1:N)
    r = rand(rng)
    if r < b/(b+d)
        N, nextID = celldivision!(alivecells, randcellidx, N, t, nextID, μ, mutationdist, rng)
    else 
        N = celldeath!(alivecells, randcellidx, N, t, rng)
    end
    return N, nextID
end

function celldivision!(alivecells::Vector{BinaryNode{SimpleCell}}, parentcellidx, N, t, nextID, μ, mutationdist, rng)

    parentcellnode = alivecells[parentcellidx] #get parent cell node
    deleteat!(alivecells, parentcellidx) #delete parent cell node from alivecells list

    #if mutations are time dependent assign mutations to parent cell and give new cells
    #no initial mutations
    childcellmuts1, childcellmuts2 = 0, 0
    if mutationdist == :fixedtimedep || mutationdist == :poissontimedep    
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

# #currently for neutral evolution only, add selection later
# function branchingupdate!(alivecells::Vector{BinaryNode{SimpleCell}}, b, d, Nmax, N, t, 
#     nextID, μ, mutationdist, rng)

#     t += exptime(rng, N * (b + d))
#     randcellidx = rand(rng, 1:N)
#     r = rand(rng)
#     if r < b/(b+d)
#         N, nextID = celldivision!(alivecells, randcellidx, N, t, nextID)
#     else 
#         N = celldeath!(alivecells, randcellidx, N, t, rng)
#     end
#     return N, t, nextID
# end

# function celldivision!(alivecells::Vector{BinaryNode{SimpleCell}}, parentcellidx, N, t, nextID)

#     parentcellnode = alivecells[parentcellidx]
#     deleteat!(alivecells, parentcellidx)
#     childcell1 = SimpleCell(nextID, true, t, 0, parentcellnode.data.clonetype)    
#     childcell2 = SimpleCell(nextID + 1, true, t, 0, parentcellnode.data.clonetype)
#     push!(alivecells, leftchild(childcell1, parentcellnode))
#     push!(alivecells, rightchild(childcell2, parentcellnode)) 

#     return N + 1, nextID + 2

# end

function celldeath!(alivecells::Vector{BinaryNode{SimpleCell}}, deadcellidx, N, t, rng::AbstractRNG)

    deadcellnode = alivecells[deadcellidx]
    deleteat!(alivecells, deadcellidx)
    leftchild(SimpleCell(deadcellnode.data.id, false, t, 0, deadcellnode.data.clonetype), deadcellnode)
    return N - 1   
end

function initialize_branching_tree(clonalmutations, initialID)
    initialcell = SimpleCell(1, true, 0.0, clonalmutations, 1)
    root = BinaryNode{SimpleCell}(initialcell)
    alivecells = [root]
    return alivecells, root, 1, initialID + 1
end

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

function deathtime(cellnode::BinaryNode)
    if AbstractTrees.has_children(cellnode)
        return cellnode.left.data.birthtime
    else
        return nothing
    end
end 

function celllifetime(cellnode::BinaryNode, tmax=nothing)
    if AbstractTrees.has_children(cellnode)
        return cellnode.left.data.birthtime - cellnode.data.birthtime
    else
        tmax = isnothing(tmax) ? age(root(cellnode)) : tmax
        return tmax - cellnode.data.birthtime

    end
end

function celllifetimes(tree; excludeliving=true)
    lifetimes = Float64[]
    if excludeliving
        for cellnode in PreOrderDFS(tree)
            if AbstractTrees.has_children(cellnode)
                push!(lifetimes, cellnode.left.data.birthtime - cellnode.data.birthtime)
            end
        end
    else
        popage = age(tree)
        for cellnode in PreOrderDFS(tree)
            push!(lifetimes, celllifetime(cellnode, popage))
        end
    end
    return lifetimes
end

function age(root::BinaryNode)
    age = 0
    for cellnode in Leaves(root)
        if cellnode.data.birthtime > age
            age = cellnode.data.birthtime
        end
    end
    return age
end

function asymmetrictree(n, b, μ, mutationdist, clonalmutations, rng::AbstractRNG=Random.GLOBAL_RNG)
    alivecells, root, N, nextID = initialize_branching_tree(clonalmutations, 1)
    t = 0.0
    for cellnode in alivecells
        t += exptime(rng, N * b)
        if N >= 2^n
            break
        else
            parentcellidx = findall(x->x==cellnode, alivecells)[1]
            N, nextID = celldivision!(alivecells, parentcellidx, N, t, nextID)
        end
    end
    addmutations!(root, μ, mutationdist, t, rng)
    return alivecells, root
end

function symmetrictree(n, b, μ, mutationdist, clonalmutations, rng::AbstractRNG=Random.GLOBAL_RNG)
    alivecells, root, N, nextID = initialize_branching_tree(clonalmutations, 1)
    t = 0.0
    for i in 1:n
        for cellnode in Leaves(root)
            t += exptime(rng, N * b)
            parentcellidx = findall(x->x==cellnode, alivecells)[1]
            N, nextID = celldivision!(alivecells, parentcellidx, N, t, nextID)
        end
    end
    addmutations!(root, μ, mutationdist, t, rng)
    return alivecells, root
end