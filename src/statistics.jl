
"""
    mutations_per_cell(population)

Calculate the number of mutations per cell in each module.
"""
mutations_per_cell(population::MultiSimulation) = map(mutations_per_cell, population)

"""
    mutations_per_cell(simulation::Simulation)
Calculate the number of mutations per cell
"""
mutations_per_cell(simulation::Simulation) = mutations_per_cell(simulation.output)

"""
    mutations_per_cell(cellmodule::CellModule)
"""
mutations_per_cell(cellmodule::CellModule) = mutations_per_cell(cellmodule.cells)

mutations_per_cell(cells::CellVector) = map(cell -> length(cell.mutations), cells)

function mutations_per_cell(root::BinaryNode{T}; includeclonal=false) where T <: AbstractTreeCell
    mutspercell = Int64[]
    for cellnode in Leaves(root)
        if isalive(cellnode.data)
            mutations = cellnode.data.mutations
            while true
                if !AbstractTrees.isroot(cellnode) && (cellnode != root|| includeclonal)
                    cellnode = cellnode.parent
                    mutations += cellnode.data.mutations
                else
                    break
                end
            end
            push!(mutspercell, mutations)
        end
    end
    return mutspercell
end

function mutations_per_cell(treemodule::TreeModule)
    mutspercell = Int64[]
    for cellnode in treemodule.cells
        mutations = cellnode.data.mutations
        while true
            if isnothing(cellnode.parent)
                break
            else
                cellnode = cellnode.parent
                mutations += cellnode.data.mutations
            end
        end
        push!(mutspercell, mutations)
    end
    return mutspercell
end

function acquired_mutations(root)
    muts = Int64[]
    for cellnode in PreOrderDFS(root)
        push!(muts, cellnode.data.mutations)
    end
    return muts
end

mutation_ids_by_cell(population::SinglelevelPopulation, idx=nothing) = mutation_ids_by_cell(population.singlemodule, idx)
mutation_ids_by_cell(cellmodule::CellModule, idx=nothing) = mutation_ids_by_cell(cellmodule.cells, idx)

function mutation_ids_by_cell(cells::Vector{Cell}, idx=nothing)
    if isnothing(idx) 
        return map(cell -> cell.mutations, cells)
    else
        return map(cell -> cell.mutations, cells[idx])
    end
end
"""
    average_mutations_per_module(population)

calculate the mean number of mutations in each module 
"""
function average_mutations_per_module(population)
    return map(average_mutations, population)
end

"""
    average_mutations(population, var=false)

calculate the mean number of mutations in whole population (and variance if variance is true)
"""
function average_mutations(population, variance=false)
    mutations = (length(cell.mutations) 
        for cellmodule in population for cell in cellmodule.cells)
    if variance
        return mean(mutations), var(mutations)
    else
        return mean(mutations)
    end
end

"""
    average_mutations(simulation::Simulation)
"""
average_mutations(simulation::Simulation) = average_mutations(simulation.output)
"""
    average_mutations(cellmodule::CellModule)
"""
function average_mutations(cellmodule::CellModule)
    return mean(length(cell.mutations) for cell in cellmodule.cells)
end

function average_mutations(treemodule::TreeModule)
    return mean(mutations_per_cell(treemodule))
end

"""
    clonal_mutations(simulation::Simulation)

Returns the number of clonal mutations, i.e. mutations present in every cell
"""
function clonal_mutations(simulation::Simulation)
    return clonal_mutations(simulation.output)
end

function clonal_mutations(population::Population)
    return map(clonal_mutations, population)
end

function clonal_mutations(population::SinglelevelPopulation)
    return clonal_mutations(population.singlemodule)
end

"""
    clonal_mutations(cellmodule::CellModule)
"""
function clonal_mutations(cellmodule::CellModule)
    return length(clonal_mutation_ids(cellmodule))
end

"""
    clonal_mutations(treemodule::TreeModule)
"""
function clonal_mutations(cellmodule::TreeModule)
    MRCA = findMRCA(cellmodule)
    return all_cell_mutations(MRCA)
end

function all_cell_mutations(cellnode::BinaryNode)
    mutations = cellnode.data.mutations
    while !isnothing(cellnode.parent)
        cellnode = cellnode.parent
        mutations += cellnode.data.mutations
    end
    return mutations
end


"""
    clonal_mutation_ids(population, idx=nothing)
"""
function clonal_mutation_ids(population, idx=nothing)
    if isnothing(idx)
        return map(clonal_mutation_ids, population)
    else
        return map(clonal_mutation_ids, population[idx])
    end
end

"""
    clonal_mutation_ids(cellmodule::CellModule)
"""
function clonal_mutation_ids(cellmodule::CellModule)
    return intersect([cell.mutations for cell in cellmodule.cells]...)
end

"""
    pairwise_fixed_differences(simulation::Simulation[, idx])
    pairwise_fixed_differences(population:Population[, idx])

Calculate the number of pairwise fixed differences between every pair of cells and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed 
cells.
"""
function pairwise_fixed_differences(simulation::Simulation, idx=nothing)
    return pairwise_fixed_differences(simulation.output, idx)
end

function pairwise_fixed_differences(population::Population{T}, idx=nothing) where T <: CellModule
    clonalmutids = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences(clonalmutids)
end

function pairwise_fixed_differences(muts::Vector{Vector{T}}) where T <: Integer
    n = length(muts)
    pfd_vec = Int64[]
    for i in 1:n
        for j in i+1:n 
            push!(pfd_vec, length(symdiff(muts[i], muts[j])))
        end
    end
    return countmap(pfd_vec)
end

function pairwise_fixed_differences(population::Population{T}, idx=nothing) where T <: TreeModule
    pfdvec = _pairwise_fixed_differences(population, idx)
    return countmap(pfdvec)
end

function _pairwise_fixed_differences(population::Population{T}, idx=nothing) where T <: TreeModule

    pfd_vec = Int64[]
    MRCA_vec = isnothing(idx) ? map(findMRCA, population) : map(findMRCA, population[idx])
    n = length(MRCA_vec)
    for i in 1:n
        for j in i+1:n
            push!(pfd_vec, pairwisedistance(MRCA_vec[i], MRCA_vec[j]))
        end
    end
    return pfd_vec
end


"""
    pairwise_fixed_differences_clonal(simulation::Simulation[, idx])
    pairwise_fixed_differences_clonal(population::Population[, idx])

Calculate the number of pairwise fixed differences between every pair of modules and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed 
modules. Also return a dict giving the number of clonal mutations
and frequency.  
"""
function pairwise_fixed_differences_clonal(simulation::Simulation, idx=nothing)
    return pairwise_fixed_differences_clonal(simulation.output, idx)
end

function pairwise_fixed_differences_clonal(population::Population{T}, idx=nothing) where T <: CellModule
    clonalmutids = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences(clonalmutids), countmap(map(length, clonalmutids))
end

function pairwise_fixed_differences_clonal(population::Population{T}, idx=nothing) where T <: TreeModule
    pfdvec, clonalmutsvec = _pairwise_fixed_differences_clonal(population, idx)
    return countmap(pfdvec), countmap(clonalmutsvec)
end

function _pairwise_fixed_differences_clonal(population::Population{T}, idx=nothing) where T <: TreeModule

    pfd_vec = Int64[]
    MRCA_vec = isnothing(idx) ? map(findMRCA, population) : map(findMRCA, population[idx])
    clonalmuts = Int64[]
    n = length(MRCA_vec)
    for i in 1:n
        push!(clonalmuts, all_cell_mutations(MRCA_vec[i]))
        for j in i+1:n
            push!(pfd_vec, pairwisedistance(MRCA_vec[i], MRCA_vec[j]))
        end
    end
    return pfd_vec, clonalmuts
end

function pairwise_differences(simulation::Simulation; idx=nothing)
    return pairwise_differences(simulation.output)
end

function pairwise_differences(population::SinglelevelPopulation{T}, idx=nothing) where T <: CellModule
    cells = population.singlemodule.cells
    if !isnothing(idx)
        cells = cells[idx]
    end
    mutation_vector = Vector{Int64}[cell.mutations for cell in cells]
    return pairwise_fixed_differences(mutation_vector)
end

function pairwise_differences(population::SinglelevelPopulation{T}, idx=nothing) where T <: TreeModule
    cells = population.singlemodule.cells
    if !isnothing(idx)
        cells = cells[idx]
    end    
    n = length(cells)
    pfd_vec = Int64[]
    for i in 1:n
        for j in i+1:n 
            push!(pfd_vec, pairwisedistance(cells[i], cells[j]))
        end
    end
    return countmap(pfd_vec)
end

function pairwisedistances(root::BinaryNode{T}, idx=nothing) where T <: AbstractTreeCell
    pfd = Int64[]
    alivecells = getalivecells(root)
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(pfd, pairwisedistance(cellnode1, cellnode2))
        end
    end
    return countmap(pfd)
end

function pairwisedistance(cellnode1::BinaryNode, cellnode2::BinaryNode)
    cellnode1 == cellnode2 && return 0
    distance = 0
    while true
        #ensure that cellnode1 is higher (or equal) level in tree than cellnode2
        if  cellnode1.data.id > cellnode2.data.id
            cellnode1, cellnode2 = cellnode2, cellnode1
        end
        #if cell nodes are siblings or both are roots add their mutations to distance and return
        if (cellnode1.parent == cellnode2.parent) || (isnothing(cellnode1.parent) && isnothing(cellnode2.parent))
            return distance + cellnode1.data.mutations + cellnode2.data.mutations
        #if cellnode1 is parent of cellnode2 add mutations of cellnode2 and return
        elseif cellnode1 == cellnode2.parent
            return distance + cellnode2.data.mutations
        #if none of these are satisfied add cellnode2 mutations and reset cellnode2 to its parent
        else
            distance += cellnode2.data.mutations
            cellnode2 = cellnode2.parent
        end
    end
end

function pairwisedistance(cell1::Cell, cell2::Cell)
    return length(symdiff(cell1.mutations, cell2.mutations))
end

"""
    pairwise_fixed_differences_matrix(population[, idx], diagonals=false)

Calculate the number of pairwise fixed differences between modules. Return an n x n matrix, 
such that the value at (i,j) is the pairwise fixed differences between modules i and j, and
there are n total modules.

If idx is given only compare specified modules, otherwise compare all modules in the 
population. If diagonals is true include comparison with self (i.e. number of fixed clonal 
mutations in each module).
"""

function pairwise_fixed_differences_matrix(simulation::Simulation, idx=nothing; diagonals=false)
    return pairwise_fixed_differences_matrix(simulation.output, idx; diagonals)
end

function pairwise_fixed_differences_matrix(population::Population{T}, idx=nothing; diagonals=false) where T<:CellModule
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_matrix(clonalmuts; diagonals)
end

# function pairwise_fixed_differences_matrix(population::SinglelevelPopulation{T}, idx=nothing; diagonals=false)
#     muts = mutation_ids_by_cell(population, idx)
#     return pairwise_fixed_differences_matrix(muts, diagonals=diagonals)
# end

function pairwise_fixed_differences_matrix(muts::Vector{Vector{T}}; diagonals=false) where T <: Integer
    n = length(muts)
    pfd = zeros(Int64, n, n)
    for i in 1:n
        if diagonals pfd[i,i] = length(muts[i]) end
        for j in i+1:n
            pfd[j, i] = length(symdiff(muts[i], muts[j]))
        end
    end
    return pfd
end

function pairwise_fixed_differences_matrix(population::Population{T}, idx=nothing) where T <: TreeModule
    MRCA_vec = isnothing(idx) ? map(findMRCA, population) : map(findMRCA, population[idx])
    n = length(MRCA_vec)
    pfd = zeros(Int64, n, n)

    for i in 1:n
        if diagonals pfd[i,i] = all_cell_mutations(MRCA_vec) end
        for j in i+1:n
            pfd[j, i] = pairwisedistance(MRCA_vec[i], MRCA_vec[j])
        end
    end
    return pfd
end


"""
    pairwise_fixed_differences_statistics(simulation::Simulation, idx=nothing)

Calculate the mean and variance of the number of pairwise fixed differences between modules. 
If idx is given only compare specified modules, otherwise compare all modules in the 
population. 
"""
function pairwise_fixed_differences_statistics(simulation::Simulation, idx=nothing)
    pairwise_fixed_differences_statistics(simulation.output, idx)
end

function pairwise_fixed_differences_statistics(simulation::Simulation, samplesize::Integer, rng)
    pairwise_fixed_differences_statistics(simulation.output, samplesize, rng)
end

function pairwise_fixed_differences_statistics(population, samplesize::Int64, rng)
    idx = if isnothing(samplesize) 
            nothing
        else
            sample(rng, 1:lastindex(population), samplesize, replace=false) 
        end
    return pairwise_fixed_differences_statistics(population, idx)
end

function pairwise_fixed_differences_statistics(population::Population{T}, idx=nothing) where T <: CellModule
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_statistics(clonalmuts)
end

function pairwise_fixed_differences_statistics(clonalmuts::Vector{Vector{T}}) where T <: Integer
    n = length(clonalmuts)
    pfd = Int64[]
    for i in 1:n
        for j in i+1:n
            push!(pfd, length(symdiff(clonalmuts[i], clonalmuts[j])))
        end
    end
    nclonalmuts = (length(cmut) for cmut in clonalmuts)
    return mean(pfd), var(pfd), mean(nclonalmuts), var(nclonalmuts)
end

function pairwise_fixed_differences_statistics(population::Population{TreeModule}, idx=nothing)
    pfd, clonalmuts = _pairwise_fixed_differences_clonal(population, idx)
    return mean(pfd), var(pfd), mean(clonalmuts), var(clonalmuts)
end

"""
    shared_fixed_mutations(population[, idx])
"""
function shared_fixed_mutations(population, idx=nothing)
    return countmap(filter!(x -> x > 0, getfixedallelefreq(population, idx)))
end

function number_modules_with_mutation(clonalmuts_by_module, mutationid)
    n = 0
    for muts in clonalmuts_by_module
        if mutationid in muts
            n += 1
        end
    end
    return n
end



"""
    time_to_MRCA(cellnode1, cellnode2, t)
Computes the time that has passed between the division time of the MRCA of the two cells
and time t.
"""
function time_to_MRCA(cellnode1, cellnode2, t)

    #ensure cellnode1 is the oldest cell
    if cellnode1.data.birthtime > cellnode2.data.birthtime
        cellnode1, cellnode2 = cellnode2, cellnode1
    end

    #if either cell is the root, it is the MRCA
    isdefined(cellnode1, :parent) || return t - endtime(cellnode1)
    isdefined(cellnode2, :parent) || return t - endtime(cellnode2)

    #if cells have the same parent, that is the MRCA
    if cellnode1.parent == cellnode2.parent
        return t - cellnode1.data.birthtime
    else
         return time_to_MRCA(cellnode1, cellnode2.parent, t)
    end
end

"""
    coalescence_times(root, [idx]; t=nothing)

Computes the coalescence time (time to MRCA) for every pair of alive cells in the root and
returns as a vector.
"""
function coalescence_times(root, idx=nothing; t=nothing)
    t = isnothing(t) ? age(root) : t
    coaltimes = Float64[]
    alivecells = getalivecells(root)
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(coaltimes,time_to_MRCA(cellnode1, cellnode2, t))
        end
    end
    return coaltimes
end

"""
    conditionalfixation_times(treemodule::TreeModule, tmin=0.0)

Compute the time to fixation for all fixed mutations in `treemodule`. Only include mutations 
    that arise after `tmin`.
"""

function conditionalfixation_times(population::Population{TreeModule{TreeCell}}, tmin=0.0)
    treeroot = getroot(reduce(
        vcat, 
        [treemodule.alivecells for treemodule in population]
    ))
    conditionalfixtimes = Float64[]
    for node in AbstractTrees.PreOrderDFS
        if node.data.birthtime >= tmin && node.data.alive
            for treemodule in population
                push!(
                    conditionalfixtimes,
                    get_conditionalfixation_times(
                        node,
                        treemodule
                    )
                )
            end
        end
    end
end

getsubclonesizes(subclones::Vector{Subclone}) = map(x -> length(x), subclones)
getsubclonesizes(population::AbstractPopulation) = getsubclonesizes(population.subclones)
getsubclonesizes(simulation::Simulation) = getsubclonesizes(simulation.output)

function getsubclonesizes(mod::AbstractModule, nsubclones)
    subclonesizes = zeros(Int64, nsubclones)
    for cell in mod.cells
        subclonesizes[getclonetype(cell)] += 1
    end
    return subclonesizes
end