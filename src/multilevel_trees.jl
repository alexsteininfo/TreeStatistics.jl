
function multilevel_simulation(::Type{T}, input::MultilevelBranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG) where T <: AbstractTreeCell

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
                rng
            )
        if length(populationtracker) != 0
            return MultiSimulation(input, populationtracker)
        end
    end
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