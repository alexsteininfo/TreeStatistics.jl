function multilevel_simulation(::Type{T}, input::S, rng::AbstractRNG=Random.GLOBAL_RNG) where {T <: AbstractTreeCell, S <: MultilevelInput}

    #if the population dies out we start a new simulation
    while true 
        populationtracker = initialize_population(
            T,
            input,
            rng
        )
        populationtracker = 
            simulate!(
                populationtracker, 
                input.maxtime, 
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
                moduleupdate = S == MultilevelMoranInput ? :moran : :branching
            )
        if length(populationtracker) != 0
            return MultiSimulation(input, populationtracker)
        end
    end
end

function multilevel_simulation_timeseries(::Type{T}, input::MultilevelMoranInput, timesteps, 
    func, rng::AbstractRNG=Random.GLOBAL_RNG) where T<:AbstractTreeCell

    populationtracker = initialize_population(
        T,
        input,
        rng
    )
    data = map(timesteps) do t
        populationtracker = simulate!(
            populationtracker, 
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
            moduleupdate=:moran

        )
        return func(populationtracker)
    end
    return data 
end

function initialize_population(::Type{T}, input, rng) where T <: AbstractTreeCell
    initialmodule = TreeModule(
        Int64[1],
        Float64[0.0],
        initialize(T, input, rng),
        CloneTracker[],
        1,
        0
    )
    return TreeModule[initialmodule]
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