function runsimulation(::Type{T}, input::S, rng::AbstractRNG=Random.GLOBAL_RNG) where {T <: AbstractTreeCell, S <: MultilevelInput}

    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            T,
            input,
            rng
        )
        nextID, nextmoduleID = 2, 2
        population, = 
            simulate!(
                population, 
                input.tmax, 
                input.maxmodules, 
                input.birthrate, 
                input.deathrate, 
                input.moranrate, 
                input.asymmetricrate,
                input.branchrate, 
                input.modulesize, 
                input.branchinitsize, 
                input.modulebranching,
                input.μ,
                input.mutationdist,
                input.moranincludeself,
                nextID,
                nextmoduleID,
                rng,
                moduleupdate = S == MultilevelBranchingMoranInput ? :moran : :branching
            )
        if length(population) != 0
            return MultiSimulation(input, population)
        end
    end
end

function runsimulation_timeseries_returnfinalpop(::Type{T}, input::S, timesteps, 
    func, rng::AbstractRNG=Random.GLOBAL_RNG) where {T <: AbstractTreeCell, S <: MultilevelInput}

    moduleupdate = S == MultilevelBranchingMoranInput ? :moran : :branching
    population = initialize_population(T, input, rng)
    nextID, nextmoduleID = 2, 2
    data = []
    t0 = 0.0
    for t in timesteps
        population, nextID, nextmoduleID = simulate!(
            population, 
            t,
            input.maxmodules, 
            input.birthrate, 
            input.deathrate, 
            input.moranrate, 
            input.asymmetricrate,
            input.branchrate, 
            input.modulesize, 
            input.branchinitsize, 
            input.modulebranching,
            input.μ,
            input.mutationdist,
            input.moranincludeself,
            nextID,
            nextmoduleID,
            rng;
            moduleupdate,
            t0
        )
        #stop branching process simulations if maximum population size is exceeded
        (moduleupdate == :branching && length(population) >= input.maxmodules) && break
        push!(data, func(population))
        t0 = t
    end
    return data, population
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