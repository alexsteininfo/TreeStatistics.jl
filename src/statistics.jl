"""
    mutations_per_cell(population)

calculate the number of mutations per cell in each module 
"""
mutations_per_cell(population::MultiSimulation) = map(mutations_per_cell, population)

"""
    mutations_per_cell(simulation::Simulation)

calculate the number of mutations per cell
"""
mutations_per_cell(simulation::Simulation) = mutations_per_cell(simulation.output)

"""
    mutations_per_cell(moduletracker::ModuleTracker)
"""
mutations_per_cell(moduletracker::ModuleTracker) = mutations_per_cell(moduletracker.cells)

mutations_per_cell(cells::Array{Cell, 1}) = map(cell -> length(cell.mutations), cells)

function mutations_per_cell(tree::BinaryNode{SimpleCell}; includeclonal=false)
    mutspercell = Int64[]
    for cellnode in Leaves(tree)
        if cellnode.data.alive
            mutations = cellnode.data.mutations
            while true
                if !AbstractTrees.isroot(cellnode) && (cellnode != tree|| includeclonal)
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

function acquired_mutations(tree)
    muts = Int64[]
    for cellnode in PreOrderDFS(tree)
        push!(muts, cellnode.data.mutations)
    end
    return muts
end

mutation_ids_by_cell(moduletracker::ModuleTracker, idx=nothing) = mutation_ids_by_cell(moduletracker.cells, idx)

function mutation_ids_by_cell(cells::Array{Cell, 1}, idx=nothing)
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
        for moduletracker in population for cell in moduletracker.cells)
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
    average_mutations(moduletracker::ModuleTracker)
"""
function average_mutations(moduletracker::ModuleTracker)
    return mean(length(cell.mutations) for cell in moduletracker.cells)
end

"""
    clonal_mutations(simulation::Simulation)

returns the number of clonal mutations, i.e. mutations present in every cell
"""
function clonal_mutations(simulation::Simulation)
    return clonal_mutations(simulation.output)
end

"""
    clonal_mutations(multisim::MultiSimulation)

returns the number of clonal mutations in each module/simulated population
"""
function clonal_mutations(population::MultiSimulation)
    return map(clonal_mutations, population)
end

"""
    clonal_mutations(output::SimulationResult)
"""
function clonal_mutations(moduletracker::ModuleTracker)
    return length(clonal_mutation_ids(moduletracker))
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
    clonal_mutation_ids(moduletracker)
"""
function clonal_mutation_ids(moduletracker::ModuleTracker)
    return intersect([cell.mutations for cell in moduletracker.cells]...)
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
    pairwise_fixed_differences(moduletracker::ModuleTracker[, idx])

See pairwise_fixed_differences(simulation::Simulation)
"""

function pairwise_fixed_differences(moduletracker::ModuleTracker, idx=nothing)
    muts = mutation_ids_by_cell(moduletracker, idx)
    return pairwise_fixed_differences(muts)
end
"""
    pairwise_fixed_differences(population[, idx])

Calculate the number of pairwise fixed differences between every pair of modules and return
as a dictionary (number differneces => frequency). If `idx` is given, only include listed 
modules.
"""

function pairwise_fixed_differences(population, idx=nothing)
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences(clonalmuts)
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

# function pairwise_fixed_differences_vector(tree::BinaryNode{SimpleCell}, pfd::Vector{Int64}=Int64[])
#     if isdefined(tree, :left) && isdefined(tree, :right)
#         for m1 in mutations_per_cell(tree.left)
#             for m2 in mutations_per_cell(tree.right)
#                 pfd = push!(pfd, m1+m2)
#             end
#         end
#     end
#     if isdefined(tree, :left)
#         pairwise_fixed_differences_vector(tree.left, pfd)
#     end
#     if isdefined(tree, :right)
#         pairwise_fixed_differences_vector(tree.right, pfd)
#     end
#     return pfd
# end
# SLOWER IMPLEMENTATION OF BELOW
# function pairwise_fixed_differences(tree::BinaryNode{SimpleCell}, idx=nothing)
#     return countmap(pairwise_fixed_differences_vector(tree))
# end

function pairwise_fixed_differences(tree::BinaryNode{SimpleCell}, idx=nothing)
    pfd = Int64[]
    alivecells = [cellnode for cellnode in Leaves(tree) if cellnode.data.alive]
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(pfd, pairwisedistance(cellnode1, cellnode2))
        end
    end
    return countmap(pfd)
end

function pairwisedistance(cellnode1::BinaryNode, cellnode2::BinaryNode, distance=0)
    if cellnode1.data.id > cellnode2.data.id
        cellnode1, cellnode2 = cellnode2, cellnode1
    end
    if cellnode1 == cellnode2
        return distance + 0
    elseif cellnode1.parent == cellnode2.parent
        return distance + cellnode1.data.mutations + cellnode2.data.mutations
    elseif isdefined(cellnode2, :parent)
         return distance + cellnode2.data.mutations + pairwisedistance(cellnode1, cellnode2.parent)
    end
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

function pairwise_fixed_differences_matrix(moduletracker::ModuleTracker, idx=nothing; diagonals=false)
    muts = mutation_ids_by_cell(moduletracker, idx)
    return pairwise_fixed_differences_matrix(muts, diagonals=diagonals)
end

function pairwise_fixed_differences_matrix(muts::Vector{Vector{Int64}}; diagonals=false)
    n = length(muts)
    pfd = zeros(Int64, n, n)
    for i in 1:n
        if diagonals pfd[i,i] = length(muts[i]) end
        for j in i+1:n
            pfd[j,i] = length(symdiff(muts[i], muts[j]))
        end
    end
    return pfd
end

function pairwise_fixed_differences_matrix(tree::BinaryNode{SimpleCell})
    alivecells = [cellnode for cellnode in Leaves(tree) if cellnode.data.alive]
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
    pairwise_fixed_differences_statistics(population[, idx], clonal=false)

Calculate the mean and variance of the number of pairwise fixed differences between modules. 
If idx is given only compare specified modules, otherwise compare all modules in the 
population. If clonal is true also calculate mean and variance of number of clonal mutations
in each module.
"""
function pairwise_fixed_differences_statistics(population, idx=nothing; clonal=true)
    clonalmuts = clonal_mutation_ids(population, idx)
    return pairwise_fixed_differences_statistics(clonalmuts; clonal=clonal)
end

function pairwise_fixed_differences_statistics(clonalmuts::Vector{Vector{Int64}}; clonal=true)
    n = length(clonalmuts)
    pfd = Int64[]
    for i in 1:n
        for j in i+1:n
            push!(pfd, length(symdiff(clonalmuts[i], clonalmuts[j])))
        end
    end
    if clonal
        nclonalmuts = (length(cmut) for cmut in clonalmuts)
        return mean(pfd), var(pfd), mean(nclonalmuts), var(nclonalmuts)
    else
        return mean(pfd), var(pfd)
    end
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
newmoduletimes(population) = sort([moduletracker.tvec[1] for moduletracker in population])

"""
    numbermodules(population, tstep)

Return number of modules in the population at times in 1:tstep:tend
"""
function numbermodules(population, tstep, tend=nothing)
    if isnothing(tend)
        tend = maximum(moduletracker.tvec[end] for moduletracker in population)
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
    tend = maximum(moduletracker.tvec[end] for moduletracker in population)
    popvec = Int64[]
    for time in 0:tstep:tend
        pop = 0
        for moduletracker in population
            N0 = 0
            for (N, t) in zip(moduletracker.Nvec, moduletracker.tvec)
                if t > time 
                    pop += N0
                    break
                elseif t == moduletracker.tvec[end]
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
    tend = maximum(moduletracker.tvec[end] for moduletracker in multisimulation)
    popvec = Float64[]
    for time in 0:tstep:tend
        pop = 0
        modules = 0
        for moduletracker in multisimulation
            if moduletracker.tvec[1] <= time 
                modules += 1
                N0 = 0 
                for (N, t) in zip(moduletracker.Nvec, moduletracker.tvec)
                    if t > time 
                        pop += N0
                        break
                    elseif t == moduletracker.tvec[end]
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
    isdefined(cellnode1, :parent) || return t - deathtime(cellnode1)
    isdefined(cellnode2, :parent) || return t - deathtime(cellnode2)

    #if cells have the same parent, that is the MRCA
    if cellnode1.parent == cellnode2.parent
        return t - cellnode1.data.birthtime
    else
         return time_to_MRCA(cellnode1, cellnode2.parent, t)
    end
end

"""
    coalescence_times(tree, [idx]; t=nothing)

Computes the coalescence time (time to MRCA) for every pair of alive cells in the tree and
returns as a vector.
"""
function coalescence_times(tree, idx=nothing; t=nothing)
    t = isnothing(t) ? age(tree) : t
    coaltimes = Float64[]
    alivecells = [cellnode for cellnode in Leaves(tree) if cellnode.data.alive]
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(coaltimes,time_to_MRCA(cellnode1, cellnode2, t))
        end
    end
    return coaltimes
end
