
"""
    true_clonal_mutations(simulation::Simulation)

returns the number of clonal mutations, i.e. mutations present in every cell
"""
function true_clonal_mutations(simulation::Simulation)
    return true_clonal_mutations(simulation.output, simulation.input.ploidy)
end

"""
    true_clonal_mutations(multisim::MultiSimulation)

returns the number of clonal mutations in each module/simulated population
"""
function true_clonal_mutations(multisim::MultiSimulation)
    ploidy = multisim.input.ploidy
    return map(output -> true_clonal_mutations(output, ploidy), multisim.output)
end

"""
    true_clonal_mutations(output::SimulationResult)
"""
function true_clonal_mutations(output::SimulationResult, ploidy=2)
    return sum(output.trueVAF .≈ 1/ploidy)
end


"""
    clonal_mutation_ids(multisim::MultiSimulation, idx=nothing)
"""
function clonal_mutation_ids(multisim::MultiSimulation, idx=nothing)
    if isnothing(idx)
        return map(clonal_mutation_ids, multisim.output)
    else
        return map(clonal_mutation_ids, multisim.output[idx])
    end
end

"""
    clonal_mutation_ids(output:SimulationResult)
"""
function clonal_mutation_ids(output::SimulationResult)
    return intersect([cell.mutations for cell in output.cells]...)
end



"""
    pairwise_fixed_differences(population::Population, idx=nothing)

returns the number of pairwise fixed differences between modules sampled from the 
population as an array. If idx is given only includes specified modules, otherwise compares
all modules in the population.
"""
function pairwise_fixed_differences(population, idx=nothing, diagonals=false)
    clonalmuts = clonal_mutation_ids(population, idx)
    n = length(clonalmuts)
    pfd = zeros(Int64, n, n)
    for i in 1:n
        if diagonals pfd[i,i] = length(clonalmuts[i]) end
        for j in i+1:n
            pfd[j,i] = length(intersect(clonalmuts[i], clonalmuts[j]))
        end
    end
    return pfd
end


"""
    newmoduletimes(population)

returns a vector of times at which a new module was formed in the population
"""
newmoduletimes(population) = newmoduletimes(population.output)
newmoduletimes(output::Vector{SimulationResult}) = sort([modout.tvec[1] for modout in output])
