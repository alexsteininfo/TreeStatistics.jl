function runsimulation(::Type{T}, input::S, rng::AbstractRNG=Random.GLOBAL_RNG) where {T <: AbstractTreeCell, S <: MultilevelInput}

    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            T,
            input,
            rng
        )
        population = 
            simulate!(
                population, 
                input.tmax, 
                input.maxmodules, 
                input.b, 
                input.d, 
                input.bdrate, 
                input.branchrate, 
                input.modulesize, 
                input.branchinitsize, 
                input.μ,
                input.mutationdist,
                rng,
                moduleupdate = S == MultilevelBranchingMoranInput ? :moran : :branching
            )
        if length(population) != 0
            return MultiSimulation(input, population)
        end
    end
end

function runsimulation_timeseries(::Type{T}, input::S, timesteps, 
    func, rng::AbstractRNG=Random.GLOBAL_RNG) where {T <: AbstractTreeCell, S <: MultilevelInput}

    population = initialize_population(
        T,
        input,
        rng
    )
    data = map(timesteps) do t
        population = simulate!(
            population, 
            t,
            input.maxmodules, 
            input.b, 
            input.d, 
            input.bdrate, 
            input.branchrate, 
            input.modulesize, 
            input.branchinitsize, 
            input.μ,
            input.mutationdist,
            rng,
            moduleupdate = S == MultilevelBranchingMoranInput ? :moran : :branching

        )
        return func(population)
    end
    return data 
end

function initialize_population(::Type{T}, input, rng) where T <: AbstractTreeCell
    return TreeModule[initialize(T, input, rng)]
end

function getnextID(population::Vector{TreeModule})
    nextID = 1
    for treemodule in population
        for cellnode in treemodule.cells
            if id(cellnode) + 1 > nextID
                nextID = id(cellnode) + 1
            end
        end
    end
    return nextID
end