"""
    mutations_per_cell(population)

calculate the number of mutations per cell in each module 
"""
mutations_per_cell(population) = map(mutations_per_cell, population)

"""
    mutations_per_cell(simulation::Simulation)

calculate the number of mutations per cell
"""
mutations_per_cell(simulation::Simulation) = mutations_per_cell(simulation.output)

"""
    mutations_per_cell(moduletracker::ModuleTracker)
"""
function mutations_per_cell(moduletracker::ModuleTracker)
    return map(cell -> length(cell.mutations), moduletracker.cells)
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
function clonal_mutations(population)
    return map(clonal_mutations, population)
end

"""
    clonal_mutations(output::SimulationResult)
"""
function clonal_mutations(moduletracker::ModuleTracker)
    return length(clonal_mutation_ids(moduletracker))
end


"""
    clonal_mutation_ids(multisim::MultiSimulation, idx=nothing)
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
    pairwise_fixed_differences(population, <idx=nothing, diagonals=flase>)

Calculate the number of pairwise fixed differences between modules. Return an n x n matrix, 
such that the value at (i,j) is the pairwise fixed differences between modules i and j, and
there are n total modules.

If idx is given only compare specified modules, otherwise compare all modules in the 
population. If diagonals is true include comparison with self (i.e. number of fixed clonal 
mutations in each module).
"""
function pairwise_fixed_differences(population; idx=nothing, diagonals=false)
    clonalmuts = clonal_mutation_ids(population, idx)
    n = length(clonalmuts)
    pfd = zeros(Int64, n, n)
    for i in 1:n
        if diagonals pfd[i,i] = length(clonalmuts[i]) end
        for j in i+1:n
            pfd[j,i] = length(symdiff(clonalmuts[i], clonalmuts[j]))
        end
    end
    return pfd
end

"""
    pairwise_fixed_differences_statistics(population, <idx=nothing, diagonals=flase>)

Calculate the mean and variance of the number of pairwise fixed differences between modules. 
If idx is given only compare specified modules, otherwise compare all modules in the 
population. If clonal is true also calculate mean and variance of number of clonal mutations
in each module.
"""
function pairwise_fixed_differences_statistics(population; idx=nothing, clonal=true)
    clonalmuts = clonal_mutation_ids(population, idx)
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
    newmoduletimes(population)

Return a vector of times at which new modules arose in the population
"""
newmoduletimes(population) = sort([moduletracker.tvec[1] for moduletracker in population])

"""
    cellpopulationsize(population, tstep)

Return number of cells in the population at times in 1:tstep:tend
"""
function cellpopulationsize(population, tstep)
    tend = maximum(moduletracker.tvec[end] for moduletracker in population)
    popvec = Int64[]
    for time in 1:tstep:tend
        pop = 0
        for moduletracker in population
            N0 = 0
            for (N, t) in zip(moduletracker.Nvec, moduletracker.tvec)
                if t > time
                    pop += N0
                    break
                else
                    N0 = N
                end
            end
        end
        push!(popvec, pop)
    end
    return collect(1:tstep:tend), popvec
end

    
