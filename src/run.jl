"""
    runsimulation(
        [::Type{T},] 
        [::Type{S},] 
        input::SimulationInput, 
        selection::AbstractSelection = NeutralSelection(),
        rng::AbstractRNG = Random.GLOBAL_RNG
    )

Run a single simulation defined by `input` and `selection`. 

### Simulation type
- If `T == Cell` a cell-type simulation is run, where mutations are represented
    by unique ids. 
- If `T == SimpleTreeCell` (default) or `T == TreeCell` a tree-type simulation is run, where 
    the number of mutations is assigned to each cell and their ancestory is tracked.

### Module structure
- If `S == WellMixed` (default) modules do not have structure
- If `S == Linear` cells are arranged in a line (not currently implemented)

"""
function runsimulation end

function runsimulation(input::SimulationInput, args...)
    return runsimulation(SimpleTreeCell, WellMixed, input, NeutralSelection(), args...) 
end

function runsimulation(input::SimulationInput, selection::AbstractSelection, args...)
    return runsimulation(SimpleTreeCell, WellMixed, input, selection, args...) 
end

function runsimulation(::Type{T}, input::SimulationInput, args...) where T <: AbstractCell
    return runsimulation(T, WellMixed, input, NeutralSelection(), args...) 
end

function runsimulation(::Type{T}, input::SimulationInput, selection::AbstractSelection, args...) where T <: AbstractCell
    return runsimulation(T, WellMixed, input, selection, args...) 
end

function runsimulation(
    ::Type{T}, 
    ::Type{S}, 
    input::SimulationInput, 
    selection=NeutralSelection()::AbstractSelection,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime, returnextinct=false
) where {T, S}

    #If T==Cell: Initially set clonalmutations = 0 and μ = 1. These are expanded later. 
    #UNLESS input.mutationdist==(:poissontimedep or :fixedtimedep) or μ <=1
    reset_mutation = reset_mutationargs(T, input)
    input_original = input
    if reset_mutation
        input = newinput(input, μ=[1], clonalmutations=0, mutationdist=[:fixed])
    end
    population = initialize_population(T, S, input; rng)
    #if the population dies out we start a new simulation (unless returnextinct=true)
    while true 
        counters = initialize_counters(population)
        population, = simulate!(population, input, selection, counters, rng; timefunc)
        if length(population) != 0 || returnextinct
            break
        else
            population = initialize_population(T, S, input; rng)
        end
    end
    #if we set μ=1 earlier expand now
    if reset_mutation
        input = input_original
        population = processresults!(population, input.μ, input.clonalmutations, rng)
    end
    return Simulation(input, population)
end

# Take a previous simulation "InitialPopulation" as input for the next simulation
# If the populatoin goes extinct, we do not start again
function runsimulation(
    InitialSimulation,
    ::Type{T}, 
    ::Type{S}, 
    input::SimulationInput, 
    selection=NeutralSelection()::AbstractSelection,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime, returnextinct=false
) where {T, S}

    #If T==Cell: Initially set clonalmutations = 0 and μ = 1. These are expanded later. 
    #UNLESS input.mutationdist==(:poissontimedep or :fixedtimedep) or μ <=1
    #reset_mutation = reset_mutationargs(T, input)
    #input_original = input
    #if reset_mutation
    #    input = newinput(input, μ=[1], clonalmutations=0, mutationdist=[:fixed])
    #end

    population = InitialSimulation.output

    counters = initialize_counters(population)
    population, = simulate!(population, input, selection, counters, rng; timefunc)
    prtinln("NewFunction")

    #if we set μ=1 earlier expand now
    #if reset_mutation
    #    input = input_original
    #    population = processresults!(population, input.μ, input.clonalmutations, rng)
    #end
    return Simulation(input, population)
end

"""
    runsimulation_timeseries_returnfinalpop(
        [::Type{T},] 
        [::Type{S},] 
        input::SimulationInput, 
        selection::AbstractSelection = NeutralSelection(),
        timesteps,
        func,
        rng::AbstractRNG = Random.GLOBAL_RNG
    )

Run a single simulation defined by `input` and `selection`. Calls `func(population)` at each
`timestep` and returns as a vector along with the final population state.

### Simulation type
- If `T == Cell` a cell-type simulation is run, where mutations are represented
    by unique ids. 
- If `T == SimpleTreeCell` (default) or `T == TreeCell` a tree-type simulation is run, where 
    the number of mutations is assigned to each cell and their ancestory is tracked.

### Module structure
- If `S == WellMixed` (default) modules do not have structure
- If `S == Linear` cells are arranged in a line (not currently implemented)
"""
function runsimulation_timeseries_returnfinalpop end

function runsimulation_timeseries_returnfinalpop(input::SimulationInput, args...)
    return runsimulation_timeseries_returnfinalpop(SimpleTreeCell, WellMixed, input, NeutralSelection(), args...) 
end

function runsimulation_timeseries_returnfinalpop(input::SimulationInput, selection::AbstractSelection, args...)
    return runsimulation_timeseries_returnfinalpop(SimpleTreeCell, WellMixed, input, selection, args...) 
end

function runsimulation_timeseries_returnfinalpop(::Type{T}, input::SimulationInput, args...) where T <: AbstractCell
    return runsimulation_timeseries_returnfinalpop(T, WellMixed, input, NeutralSelection(), args...) 
end

function runsimulation_timeseries_returnfinalpop(::Type{T}, input::SimulationInput, selection::AbstractSelection, args...) where T <: AbstractCell
    return runsimulation_timeseries_returnfinalpop(T, WellMixed, input, selection, args...) 
end

function runsimulation_timeseries_returnfinalpop(
    ::Type{T}, 
    ::Type{S}, 
    input, 
    selection, 
    timesteps, 
    func, 
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime, 
    returnextinct=false
) where {T, S}

    population = initialize_population(T, S, input; rng)
    data = []
    moduleupdate = getmoduleupdate(input)
    while true
        counters = initialize_counters(population)
        t0 = 0.0
        for t in timesteps
            simulate!(population, input, selection, counters, rng; timefunc, t0, tmax=t)
            #stop branching process simulations if maximum population size is exceeded
            if moduleupdate == :branching && length(population) >= input.maxmodules
                break
            end
            push!(data, func(population))
            t0 = t
        end
        if length(population) != 0 || returnextinct
            break
        else
            population = initialize_population(T, S, input; rng)
            data = []
        end
    end
    return data, Simulation(input, population)
end


"""
    runsimulation_timeseries(
        [::Type{T},] 
        [::Type{S},] 
        input::SimulationInput, 
        selection::AbstractSelection = NeutralSelection(),
        timesteps,
        func,
        rng::AbstractRNG = Random.GLOBAL_RNG
    )

Run a single simulation defined by `input` and `selection`. Calls `func(population)` at each
`timestep` and returns as a vector.

### Simulation type
- If `T == Cell` a cell-type simulation is run, where mutations are represented
    by unique ids. 
- If `T == SimpleTreeCell` (default) or `T == TreeCell` a tree-type simulation is run, where 
    the number of mutations is assigned to each cell and their ancestory is tracked.

### Module structure
- If `S == WellMixed` (default) modules do not have structure
- If `S == Linear` cells are arranged in a line (not currently implemented)
"""
function runsimulation_timeseries end

function runsimulation_timeseries(input::SimulationInput, args...)
    return runsimulation_timeseries(SimpleTreeCell, WellMixed, input, NeutralSelection(), args...) 
end

function runsimulation_timeseries(input::SimulationInput, selection::AbstractSelection, args...)
    return runsimulation_timeseries(SimpleTreeCell, WellMixed, input, selection, args...) 
end

function runsimulation_timeseries(::Type{T}, input::SimulationInput, args...) where T <: AbstractCell
    return runsimulation_timeseries(T, WellMixed, input, NeutralSelection(), args...) 
end

function runsimulation_timeseries(::Type{T}, input::SimulationInput, selection::AbstractSelection, args...) where T <: AbstractCell
    return runsimulation_timeseries(T, WellMixed, input, selection, args...) 
end

function runsimulation_timeseries(
    ::Type{T}, 
    ::Type{S}, 
    input, 
    selection, 
    timesteps, 
    func, 
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime, 
    returnextinct=false
) where {T, S}    
return runsimulation_timeseries_returnfinalpop(T, S, input, selection, timesteps, func, rng; 
    timefunc, returnextinct)[1] 
end

function runsimulation_condfixtime(input::MultilevelInput, rng; timefunc=exptime, returnextinct=false)
    return runsimulation_condfixtime(Cell, WellMixed, input, NeutralSelection(), rng; timefunc, returnextinct) 
end


function runsimulation_condfixtime(
    ::Type{Cell}, 
    ::Type{S}, 
    input::MultilevelInput, 
    selection::NeutralSelection, 
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime, 
    returnextinct=false
) where S <: ModuleStructure

    #if the population dies out we start a new simulation
    population = initialize_population(Cell, S, input; rng)
    condfixtimes = (homeostatic=Vector{Float64}[[]], growing=Vector{Float64}[[]])
    while true 
        counters = initialize_counters(population)
        population, condfixtimes = simulate_condfixtime!(population, input, selection, counters, rng; timefunc)
        if length(population) != 0 || returnextinct
            break
        else
            population = initialize_population(T, S, input; rng)
        end
    end
    if length(population) != 0
        return condfixtimes, Simulation(input, population)
    end
end


getinputrates(input::BranchingInput) = input.birthrate, input.deathrate, 0.0, 0.0
getinputrates(input::MoranInput) = 0.0, 0.0, input.moranrate, 0.0
getinputrates(input::BranchingMoranInput) = input.birthrate, input.deathrate, input.moranrate, 0.0

initialize_counters(population::MultilevelPopulation) = (getnextID(population), length(population) + 1)
initialize_counters(population::SinglelevelPopulation) = getnextID(population)


function reset_mutationargs(::Type{Cell}, input)
    return (!(
        (:poissontimedep in input.mutationdist) || (:fixedtimedep in input.mutationdist)
             || (sum(input.μ) <= 1)))
end

function reset_mutationargs(::Type{T}, input) where T <: AbstractTreeCell
    return false
end