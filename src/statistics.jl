
"""
    mutations_per_cell(simulation::Simulation)
    mutations_per_cell(population::AbstractPopulation)
    mutations_per_cell(module::AbstractModule)

Calculate the number of mutations per cell in each module.
"""
mutations_per_cell(simulation::Simulation) = mutations_per_cell(simulation.output)
mutations_per_cell(population::MultilevelPopulation) = map(mutations_per_cell, population)

function mutations_per_cell(population::SinglelevelPopulation)
    return mutations_per_cell(population.singlemodule)
end

mutations_per_cell(cellmodule::CellModule) = mutations_per_cell(cellmodule.cells)
mutations_per_cell(cells::CellVector) = map(cell -> length(cell.mutations), cells)

function mutations_per_cell(
    root::BinaryNode{T}; includeclonal=false
) where T <: AbstractTreeCell
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


"""
    mutation_ids_by_cell(simulation::Simulation)
    mutation_ids_by_cell(population::AbstractPopulation)
    mutation_ids_by_cell(module::AbstractModule)

Return the ids of clonal mutations in each module (Cell type simulations only).
"""
mutation_ids_by_cell(simulation::Simulation)

function mutation_ids_by_cell(population::MultilevelPopulation)
    return map(mutation_ids_by_cell, population)
end

function mutation_ids_by_cell(population::SinglelevelPopulation)
    return mutation_ids_by_cell(population.singlemodule)
end

mutation_ids_by_cell(cellmodule::CellModule) = mutation_ids_by_cell(cellmodule.cells)
mutation_ids_by_cell(cells::Vector{Cell}) = map(cell -> cell.mutations, cells)

"""
    average_mutations(simulation::Simulation)
    average_mutations(population::AbstractPopulation)
    average_mutations(module::AbstractModule)

Calculate the mean number of mutations averaged over all cells in population.
"""
average_mutations(simulation::Simulation) = average_mutations(simulation.output)

function average_mutations(population::MultilevelPopulation)
    mutations = (m
        for mod in population for m in mutations_per_cell(mod))
    return mean(mutations)
end

function average_mutations(population::SinglelevelPopulation)
    return average_mutations(population.singlemodule)
end

function average_mutations(cellmodule::CellModule)
    return mean(length(cell.mutations) for cell in cellmodule.cells)
end

function average_mutations(treemodule::TreeModule)
    return mean(mutations_per_cell(treemodule))
end

"""
    average_mutations_per_module(simulation::Simulation)
    average_mutations_per_module(population::MultilevelPopulation)
"""
function average_mutations_per_module(simulation::Simulation)
    return average_mutations_per_module(simulation.output)
end

function average_mutations_per_module(population::MultilevelPopulation)
    return map(average_mutations, population)
end

"""
    clonal_mutations(simulation::Simulation)
    clonal_mutations(population::AbstractPopulation)
    clonal_mutations(module::AbstractModule)

Return the number of clonal mutations, i.e. mutations present in every cell, for each module.
"""
function clonal_mutations(simulation::Simulation)
    return clonal_mutations(simulation.output)
end

function clonal_mutations(population::MultilevelPopulation)
    return map(clonal_mutations, population)
end

function clonal_mutations(population::SinglelevelPopulation)
    return clonal_mutations(population.singlemodule)
end

function clonal_mutations(cellmodule::CellModule)
    return length(clonal_mutation_ids(cellmodule))
end

function clonal_mutations(cellmodule::TreeModule)
    MRCA = findMRCA(cellmodule)
    return all_cell_mutations(MRCA)
end

"""
    all_cell_mutations(node::BinaryNode)

Return the total number of mutations for the cell at `node`, including those inherited.
"""
function all_cell_mutations(node::BinaryNode)
    mutations = node.data.mutations
    while !isnothing(node.parent)
        node = node.parent
        mutations += node.data.mutations
    end
    return mutations
end


"""
    clonal_mutation_ids(population::Simulation[, idx])
    clonal_mutation_ids(population::AbstractPopulation[, idx])
    clonal_mutation_ids(module::CellModule)


Return the ids of clonal mutations in each module. If `idx` is given, only include listed
modules. Only applies to `Cell` type simulations.
"""
clonal_mutation_ids(simulation, idx=nothing) = clonal_mutation_ids(simulation.output, idx)

function clonal_mutation_ids(population::MultilevelPopulation, idx=nothing)
    if isnothing(idx)
        return map(clonal_mutation_ids, population)
    else
        return map(clonal_mutation_ids, population[idx])
    end
end

function clonal_mutation_ids(population::SinglelevelInput)
    return clonal_mutation_ids(population.singlemodule)
end

function clonal_mutation_ids(cellmodule::CellModule)
    return intersect([cell.mutations for cell in cellmodule.cells]...)
end

"""
    pairwise_fixed_differences(simulation::Simulation[, idx])
    pairwise_fixed_differences(population::MultilevelPopulation[, idx])

Calculate the number of pairwise fixed differences between every pair of cells and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed
cells.
"""
function pairwise_fixed_differences(simulation::Simulation, idx=nothing)
    return pairwise_fixed_differences(simulation.output, idx)
end

function pairwise_fixed_differences(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing
) where T <: CellModule
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

function pairwise_fixed_differences(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing
) where T <: TreeModule
    pfdvec = _pairwise_fixed_differences(population, idx)
    return countmap(pfdvec)
end

function _pairwise_fixed_differences(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing
) where T <: TreeModule

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
    pairwise_fixed_differences_clonal(population::MultilevelPopulation[, idx])

Calculate the number of pairwise fixed differences between every pair of modules and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed
modules. Also return a dict giving the number of clonal mutations
and frequency.
"""
function pairwise_fixed_differences_clonal(simulation::Simulation, idx=nothing)
    return pairwise_fixed_differences_clonal(simulation.output, idx)
end

function pairwise_fixed_differences_clonal(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing
) where T <: CellModule
    clonalmutids = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences(clonalmutids), countmap(map(length, clonalmutids))
end

function pairwise_fixed_differences_clonal(
    population::Union{Population{T}, PopulationWithQuiescence{T}}, idx=nothing
) where T <: TreeModule
    pfdvec, clonalmutsvec = _pairwise_fixed_differences_clonal(population, idx)
    return countmap(pfdvec), countmap(clonalmutsvec)
end

function _pairwise_fixed_differences_clonal(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing
) where T <: TreeModule
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

"""
pairwise_differences(simulation::Simulation[, idx])
pairwise_differences(population::SinglelevelPopulation[, idx])

Calculate the number of pairwise fixed differences between every pair of cells and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed
cells.
"""
function pairwise_differences(simulation::Simulation, idx=nothing)
    return pairwise_differences(simulation.output, idx)
end

function pairwise_differences(population::SinglelevelPopulation, idx=nothing)
    return pairwise_differences(population.singlemodule.cells, idx)
end

function pairwise_differences(cells::AbstractCellVector, idx=nothing)
    if !isnothing(idx)
        cells = cells[idx]
    end
    mutation_vector = Vector{Int64}[cell.mutations for cell in cells]
    return pairwise_fixed_differences(mutation_vector)
end


function pairwise_differences(cells::AbstractTreeCellVector, idx=nothing)
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

"""
    pairwisedistance(node1::BinaryNode, node2::BinaryNode)
    pairwisedistance(cell1::Cell, cell2::Cell)

Return the number of pairwise differences between two cells.
"""
function pairwisedistance(cellnode1::BinaryNode, cellnode2::BinaryNode)
    cellnode1 == cellnode2 && return 0
    distance = 0
    while true
        #ensure that cellnode1 is higher (or equal) level in tree than cellnode2
        if  cellnode1.data.id > cellnode2.data.id
            cellnode1, cellnode2 = cellnode2, cellnode1
        end
        #if cell nodes are siblings or both are roots add their mutations to distance and return
        if (
            (cellnode1.parent == cellnode2.parent) ||
            (isnothing(cellnode1.parent) && isnothing(cellnode2.parent))
        )
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
    pairwise_fixed_differences_matrix(simulation::Simulation[, idx], diagonals=false)
    pairwise_fixed_differences_matrix(population::MultilevelPopulation[, idx],
    diagonals=false)

Calculate the number of pairwise fixed differences between modules. Return an n x n matrix,
such that the value at (i,j) is the pairwise fixed differences between modules i and j, and
there are n total modules.

If `idx` is given only compare specified modules, otherwise compare all modules in the
population. If `diagonals==true` include comparison with self (i.e. number of fixed clonal
mutations in each module).
"""

function pairwise_fixed_differences_matrix(
    simulation::Simulation,
    idx=nothing;
    diagonals=false
)
    return pairwise_fixed_differences_matrix(simulation.output, idx; diagonals)
end

function pairwise_fixed_differences_matrix(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing;
    diagonals=false
) where T<:CellModule
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_matrix(clonalmuts; diagonals)
end

function pairwise_fixed_differences_matrix(
    muts::Vector{Vector{T}};
    diagonals=false
) where T <: Integer
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

function pairwise_fixed_differences_matrix(
    population::Union{Population{T}, PopulationWithQuiescence{T}},
    idx=nothing
) where T <: TreeModule
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
    pairwise_fixed_differences_statistics(simulation::Simulation[, samplesize, idx])
    pairwise_fixed_differences_statistics(population::MultilevelPopulation[, samplesize,
        idx])

Calculate the mean and variance of the number of pairwise fixed differences between modules.
If `idx` is given only compare specified modules, otherwise compare all modules in the
population. If `samplesize` select modules at random to include (not compatible with `idx`).
"""
function pairwise_fixed_differences_statistics(simulation::Simulation, idx=nothing)
    pairwise_fixed_differences_statistics(simulation.output, idx)
end

function pairwise_fixed_differences_statistics(
    simulation::Simulation,
    samplesize::Integer,
    rng
)
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

function pairwise_fixed_differences_statistics(
    population::Population{T},
    idx=nothing
) where T <: CellModule
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_statistics(clonalmuts)
end

function pairwise_fixed_differences_statistics(
    clonalmuts::Vector{Vector{T}}
) where T <: Integer
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

function pairwise_fixed_differences_statistics(
    population::Population{TreeModule},
    idx=nothing
)
    pfd, clonalmuts = _pairwise_fixed_differences_clonal(population, idx)
    return mean(pfd), var(pfd), mean(clonalmuts), var(clonalmuts)
end

"""
    shared_fixed_mutations(simulation::Simulation[, idx])
    shared_fixed_mutations(population::MultilevelPopulation[, idx])

For every fixed mutation compute the number of modules that it is fixed in and return as
a `Dict`.
"""
function shared_fixed_mutations(simulation::Simulation, idx=nothing)
    return shared_fixed_mutations(simulation.output, idx)
end

function shared_fixed_mutations(population, idx=nothing)
    return countmap(filter!(x -> x > 0, getfixedallelefreq(population, idx)))
end

"""
    time_to_MRCA(node1::BinaryNode, node2::BinaryNode, t)

Computes the time that has passed between the division time of the MRCA of the two cells
and time t. `Tree` type simulations only
"""
function time_to_MRCA(node1, node2, t)

    #ensure node1 is the oldest cell
    if node1.data.birthtime > node2.data.birthtime
        node1, node2 = node2, node1
    end

    #if either node is the root, it is the MRCA
    isdefined(node1, :parent) || return t - endtime(node1)
    isdefined(node2, :parent) || return t - endtime(node2)

    #if nodes have the same parent, that is the MRCA
    if node1.parent == node2.parent
        return t - node1.data.birthtime
    else
         return time_to_MRCA(node1, node2.parent, t)
    end
end

"""
    coalescence_times(root, [idx]; t=nothing)

Computes the coalescence time (time to MRCA) for every pair of alive cells in the `root`.
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
