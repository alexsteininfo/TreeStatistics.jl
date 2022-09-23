
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

mutations_per_cell(cells::Vector{Cell}) = map(cell -> length(cell.mutations), cells)

function mutations_per_cell(root::BinaryNode{T}; includeclonal=false) where T <: AbstractTreeCell
    mutspercell = Int64[]
    for cellnode in Leaves(root)
        if alive(cellnode.data)
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

function acquired_mutations(root)
    muts = Int64[]
    for cellnode in PreOrderDFS(root)
        push!(muts, cellnode.data.mutations)
    end
    return muts
end

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

"""
    clonal_mutations(multisim::MultiSimulation)
"""
function clonal_mutations(population::MultiSimulation)
    return map(clonal_mutations, population)
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
    clonal_mutation_ids(cellmodule)
"""
function clonal_mutation_ids(cellmodule::CellModule)
    return intersect([cell.mutations for cell in cellmodule.cells]...)
end

"""
    pairwise_fixed_differences(simulation::Simulation[, idx])

Calculate the number of pairwise fixed differences between every pair of cells and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed 
cells.
"""

function pairwise_fixed_differences(simulation::Simulation, idx=nothing)
    return pairwise_fixed_differences(simulation.output, idx)
end

"""
    pairwise_fixed_differences(cellmodule::CellModule[, idx])

See pairwise_fixed_differences(simulation::Simulation)
"""

function pairwise_fixed_differences(cellmodule::CellModule, idx=nothing)
    muts = mutation_ids_by_cell(cellmodule, idx)
    return pairwise_fixed_differences(muts)
end
"""
    pairwise_fixed_differences(population[, idx])

Calculate the number of pairwise fixed differences between every pair of modules and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed 
modules. If `clonal` then additionally return a dict giving the number of clonal mutations
and frequency.  
"""

function pairwise_fixed_differences_clonal(population, idx=nothing)
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences(clonalmuts), countmap(map(length, clonalmuts))
end

function pairwise_fixed_differences(muts::Vector{Vector{Int64}})
    n = length(muts)
    pfd_vec = Int64[]
    for i in 1:n
        for j in i+1:n 
            push!(pfd_vec, length(symdiff(muts[i], muts[j])))
        end
    end
    return countmap(pfd_vec)
end

function pairwise_fixed_differences_clonal(population::Union{MultiSimulation{T, S}, Vector{T}}, 
    idx=nothing) where {T, S <: TreeModule}

    pfdvec, clonalmutsvec = _pairwise_fixed_differences_clonal(population, idx)
    return countmap(pfdvec), countmap(clonalmutsvec)
end

function _pairwise_fixed_differences_clonal(population::Union{MultiSimulation{T, S}, Vector{T}}, 
    idx=nothing) where {T, S <: TreeModule}

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

function pairwise_fixed_differences(module1::TreeModule, module2::TreeModule)
    return pairwisedistance(findMRCA(module1), findMRCA(module2))
end

function pairwise_fixed_differences(root::BinaryNode{T}, idx=nothing) where T <: AbstractTreeCell
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

function pairwise_fixed_differences(sampledcells::Vector{BinaryNode{TreeCell}})
    pfd = Int64[]
    while length(sampledcells) > 1
        cellnode1 = popfirst!(sampledcells)
        for cellnode2 in sampledcells
            push!(pfd, pairwisedistance(cellnode1, cellnode2))
        end
    end
    return countmap(pfd)
end


function pairwisedistance_recursive(cellnode1::BinaryNode, cellnode2::BinaryNode)
    if cellnode1.data.id > cellnode2.data.id
        cellnode1, cellnode2 = cellnode2, cellnode1
    end
    if cellnode1 == cellnode2
        return 0
    elseif cellnode1.parent == cellnode2.parent
        return cellnode1.data.mutations + cellnode2.data.mutations
    elseif cellnode1 == cellnode2.parent
        return cellnode2.data.mutations
    elseif !isnothing(cellnode2.parent)
         return cellnode2.data.mutations + pairwisedistance(cellnode1, cellnode2.parent)
    end
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
function pairwise_fixed_differences_matrix(population, idx=nothing; diagonals=false)
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_matrix(clonalmuts, diagonals=diagonals)
end

function pairwise_fixed_differences_matrix(simulation::Simulation, idx=nothing; diagonals=false)
    return pairwise_fixed_differences_matrix(simulation.output, idx, diagonals=diagonals)
end

function pairwise_fixed_differences_matrix(cellmodule::CellModule, idx=nothing; diagonals=false)
    muts = mutation_ids_by_cell(cellmodule, idx)
    return pairwise_fixed_differences_matrix(muts, diagonals=diagonals)
end

function pairwise_fixed_differences_matrix(muts::Vector{Vector{Int64}}; diagonals=false)
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

function pairwise_fixed_differences_matrix(population::MultiSimulation{T, S}, idx=nothing; 
    diagonals=false) where {T, S <: TreeModule}

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

function pairwise_fixed_differences_matrix(root::BinaryNode{T}) where T <: AbstractTreeCell
    alivecells = getalivecells(root)
    n = length(alivecells)
    pfd = zeros(Int64, n, n)
    for i in 1:n
        cellnode1 = popfirst!(alivecells)
        for (j, cellnode2) in enumerate(alivecells)
            pfd[i+j, i] = pairwisedistance(cellnode1, cellnode2)
        end
    end
    return pfd
end

"""
    pairwise_fixed_differences_statistics(population, samplesize::Int64, rng)

Calculate the mean and variance of the number of pairwise fixed differences between modules. 
If idx is given only compare specified modules, otherwise compare all modules in the 
population. If clonal is true also calculate mean and variance of number of clonal mutations
in each module.
"""
function pairwise_fixed_differences_statistics(population, samplesize::Int64, rng)
    idx = begin
        if isnothing(samplesize) 
            nothing
        else
            sample(rng, 1:lastindex(population), samplesize, replace=false) 
        end
    end
    return pairwise_fixed_differences_statistics(population, idx)
end

function pairwise_fixed_differences_statistics(population, idx=nothing)
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_statistics(clonalmuts)
end

function pairwise_fixed_differences_statistics(clonalmuts::Vector{Vector{Int64}})
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

function pairwise_fixed_differences_statistics(population::Union{MultiSimulation{T, S}, Vector{S}}, 
    idx=nothing) where {T, S <: TreeModule}

    pfd, clonalmuts = _pairwise_fixed_differences_clonal(population, idx)
    return mean(pfd), var(pfd), mean(clonalmuts), var(clonalmuts)
end

"""
    shared_fixed_mutations(population[, idx])
"""
function shared_fixed_mutations(population, idx=nothing)
    clonalmuts = clonal_mutation_ids(population, idx)
    return shared_fixed_mutations(clonalmuts)
end

function shared_fixed_mutations(clonalmuts::Vector{Vector{Int64}})
    clonalmuts_vec = reduce(union, clonalmuts)
    nclonalmuts = map(
        x -> number_modules_with_mutation(clonalmuts, x), clonalmuts_vec
    )
    return countmap(nclonalmuts)
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
    newmoduletimes(population)

Return a vector of times at which new modules arose in the population
"""
newmoduletimes(population) = sort([cellmodule.tvec[1] for cellmodule in population])

"""
    numbermodules(population, tstep)

Return number of modules in the population at times in 1:tstep:tend
"""
function numbermodules(population, tstep, tend=nothing)
    if isnothing(tend)
        tend = maximum(cellmodule.tvec[end] for cellmodule in population)
    end
    newmodtimes = newmoduletimes(population)
    times = collect(0:tstep:tend)
    nmodules = Int64[]
    for t in times
        push!(nmodules, sum(newmodtimes .<= t))
    end
    return times, nmodules

end
"""
    cellpopulationsize(population, tstep)

Return number of cells in the population at times in 1:tstep:tend
"""
function cellpopulationsize(population, tstep)
    tend = maximum(cellmodule.tvec[end] for cellmodule in population)
    popvec = Int64[]
    for time in 0:tstep:tend
        pop = 0
        for cellmodule in population
            N0 = 0
            for (N, t) in zip(cellmodule.Nvec, cellmodule.tvec)
                if t > time 
                    pop += N0
                    break
                elseif t == cellmodule.tvec[end]
                    pop += N
                    break
                else
                    N0 = N
                end
            end
        end
        push!(popvec, pop)
    end
    return collect(0:tstep:tend), popvec
end

"""
    meanmodulesize(multisimulation, tstep)

Return mean modulesize in the multisimulation at times in 1:tstep:tend
"""
function meanmodulesize(multisimulation, tstep)
    tend = maximum(cellmodule.tvec[end] for cellmodule in multisimulation)
    popvec = Float64[]
    for time in 0:tstep:tend
        pop = 0
        modules = 0
        for cellmodule in multisimulation
            if cellmodule.tvec[1] <= time 
                modules += 1
                N0 = 0 
                for (N, t) in zip(cellmodule.Nvec, cellmodule.tvec)
                    if t > time 
                        pop += N0
                        break
                    elseif t == cellmodule.tvec[end]
                        pop += N
                        break
                    else
                        N0 = N
                    end
                end
            end
        end
        push!(popvec, pop/modules)
    end
    return collect(0:tstep:tend), popvec
end

"""
    time_to_MRCA(cellnode1, cellnode2, t)
Computes the time thaat has passed between the division time of the MRCA of the two cells
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
