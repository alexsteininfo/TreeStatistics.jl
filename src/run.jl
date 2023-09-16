"""
    runsimulation(::Type{T}, ::Type{S}, input::SimulationInput, rng::AbstractRNG = MersenneTwister)

    Run a single simulation.
    Simulate a growing population of modules, until it reaches `maxmodules`.
    
    Start with a single module comprised of a single cell that grows via a branching process  
    to size `input.modulesize`, then switches to a Moran process. 
    New modules are created with rate `input.branchrate`, by sampling cells from the parent 
    module. Return output as a MultiSimulation.

    Simulation is implemented by a Gillespie algorithm and runs until the number of modules 
    exceeds input.maxmodules or time exceeds input.tmax.


"""
function runsimulation(input::SimulationInput, args...)
    return runsimulation(SimpleTreeCell, WellMixed, input, args...) 
end

function runsimulation(::Type{T}, input::SimulationInput, args...) where T <: AbstractCell
    return runsimulation(T, WellMixed, input, args...) 
end

function runsimulation(
    ::Type{T}, 
    ::Type{S}, 
    input::MultilevelInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG
) where {T, S}

    #If T==Cell: Initially set clonalmutations = 0 and μ = 1. These are expanded later. 
    #UNLESS input.mutationdist=:poissontimedep or :fixedtimedep
    μ, clonalmutations, mutationdist = getmutationargs(T, input)

    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            T, S, 
            clonalmutations, 
            getNinit(input),
            input.birthrate,
            input.deathrate,
            input.moranrate,
            input.asymmetricrate;
            rng)
        nextID, nextmoduleID = 2, 2
        population, = 
            simulate!(
                population, 
                input.tmax, 
                input.maxmodules, 
                input.branchrate, 
                input.modulesize, 
                input.branchinitsize, 
                input.modulebranching,
                μ,
                mutationdist,
                input.moranincludeself,
                nextID,
                nextmoduleID,
                rng,
                moduleupdate=getmoduleupdate(input)
            )
        if length(population) != 0
            break
        end
    end
    #if we set μ=1 earlier expand now
    if μ !== input.μ
        population = processresults!(population, input.μ, input.clonalmutations, rng)
    end

    return MultiSimulation(input, population)
end

function runsimulation(::Type{Cell}, ::Type{S}, input::SinglelevelInput, rng::AbstractRNG=Random.GLOBAL_RNG) where S <: ModuleStructure
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later. 
    #UNLESS input.mutationdist=:poissontimedep or :fixedtimedep

    μ, clonalmutations, mutationdist = input.μ, input.clonalmutations, input.mutationdist
    if !(input.mutationdist ∈ (:poissontimedep, :fixedtimedep))
        input = newinput(input, μ=1, clonalmutations=0, mutationdist=:fixed)
    end

    cellmodule = initialize(Cell, S, input.clonalmutations, getNinit(input); rng)

    simulate!(cellmodule, input, rng)
    
    #If we set μ=1 etc we need to add proper mutations now (also remove undetectable subclones).
    if !(input.mutationdist ∈ (:poissontimedep, :fixedtimedep) )  
        processresults!(cellmodule, μ, clonalmutations, rng; mutationdist)
        input = newinput(input; μ, clonalmutations, mutationdist)
    #Otherwise if mutation accumulation is time dependent, add the final mutations.
    else
        final_timedep_mutations!(cellmodule::CellModule, input.μ, input.mutationdist, rng)
    end
    return Simulation(input, cellmodule)
end

function runsimulation(::Type{T}, ::Type{S}, input::SinglelevelInput, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, returnextinct=false) where {T <: AbstractTreeCell, S <: ModuleStructure}
    while true
        treemodule = initialize(T, S, input.clonalmutations, getNinit(input); rng)
        simulate!(treemodule, input, rng; timefunc)
        if length(treemodule) > 0 || returnextinct
            return Simulation(input, treemodule)
        end
    end

end

function runsimulation_timeseries_returnfinalpop(::Type{T}, ::Type{S}, input, timesteps, func, rng::AbstractRNG=Random.GLOBAL_RNG) where {T, S}
    population = initialize_population(
        T,
        S,
        input.clonalmutations,
        getNinit(input);
        rng
    )
    data = []
    t0 = 0.0
    nextID, nextmoduleID = 2, 2
    moduleupdate = getmoduleupdate(input)
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
            moduleupdate=moduleupdate,
            t0
        )
        #stop branching process simulations if maximum population size is exceeded
        (moduleupdate == :branching && length(population) >= input.maxmodules) && break
        push!(data, func(population))
        t0 = t
    end
    return data, population
end

"""
    runsimulation_timeseries(::Type{T}, input::MultilevelInput, timesteps, func, rng::AbstractRNG=Random.GLOBAL_RNG) where T

TBW
"""
function runsimulation_timeseries(::Type{T}, ::Type{S}, input::SimulationInput, timesteps, func, rng::AbstractRNG=Random.GLOBAL_RNG) where {T, S}
    return runsimulation_timeseries_returnfinalpop(T, S, input, timesteps, func, rng)[1] 
end

function runsimulation_timeseries(::Type{T}, input::SimulationInput, timesteps, func, rng::AbstractRNG=Random.GLOBAL_RNG) where {T}
    return runsimulation_timeseries(T, WellMixed, input, timesteps, func, rng)[1] 
end

getmoduleupdate(::MultilevelBranchingInput) = :branching
getmoduleupdate(::MultilevelBranchingMoranInput) = :moran
getmoduleupdate(::MultilevelMoranInput) = :moran

function getmutationargs(::Type{Cell}, input) 
    if input.mutationdist ∈ (:poissontimedep, :fixedtimedep)
        return input.μ, input.clonalmutations, input.mutationdist
    else
        return 1, 0, :fixed
    end
end

function getmutationargs(::Type{T}, input) where T <: AbstractTreeCell
    return input.μ, input.clonalmutations, input.mutationdist
end


function run_multilevel_from_file(filename, outputdir)
    inputdict = JSON.parsefile(filename, dicttype=Dict{Symbol, Any})
    input = loadinput(MultilevelInput, inputdict[:input])
    seeds = inputdict[:seeds]
    output = inputdict[:output]
    if inputdict[:repeat] > 1
        for i in 1:inputdict[:repeat]
            rng, seed = seeds !== nothing ? (MersenneTwister(seeds[i]), seeds[i]) : (MersenneTwister(), nothing)
            population = runsimulation(input, rng, Symbol(inputdict[:simtype]))
            saveoutput(population, output, outputdir, seed, i)
        end
    else
        rng = MersenneTwister(seeds)
        population = runsimulation(input, rng, Symbol(inputdict[:simtype]))
        saveoutput(population, output, outputdir, seeds, 1)
    end
end

function run_multilevel_from_file(filename, outputdir, id)
    inputdict = JSON.parsefile(filename, dicttype=Dict{Symbol, Any})
    input = loadinput(MultilevelInput, inputdict[:input])
    seed = (inputdict[:seeds])[id]
    output = inputdict[:output]
    rng = MersenneTwister(seed)
    population = runsimulation(input, rng, Symbol(inputdict[:simtype]))
    saveoutput(population, output, outputdir, seed, id)
end