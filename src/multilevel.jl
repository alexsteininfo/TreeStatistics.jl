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

"""
    runsimulation(input::MultilevelInput, rng::AbstractRNG = MersenneTwister; 
        <keyword arguments>)

    Simulate a growing population of modules, until it reaches `maxmodules`.
    
    Start with a single module comprised of a single cell that grows via a branching process  
    to size `input.modulesize`, then switches to a Moran process. 
    New modules are created with rate `input.branchrate`, by sampling cells from the parent 
    module. Return output as a MultiSimulation.

    - if `simtype == :normal`, simulation is
    implemented by a Gillespie algorithm and runs until the number of modules exceeds 
    input.maxmodules or time exceeds input.tmax.

    - if `simtype == :fixedtime`, a list of modules is created (initially of length 1)
    and each module is simulated independently until input.tmax is reached. New modules 
    are added to the end of the list. This implementation is faster, but only works if we 
    end the simulation at a fixed time (not fixed size). If the list of modules is bigger 
    than input.maxmodules and error is thrown.

"""
function runsimulation(input::MultilevelInput, args...)
    return runsimulation(Cell, input, args...) 
end

function runsimulation(::Type{Cell}, input::MultilevelBranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG, simtype=:normal) 
    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            input.modulesize, 
            clonalmutations=0
        )
        if simtype == :normal || simtype == "normal"
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
                    1, #assume a single mutation at each division and add distribution of
                    :fixed, #mutations later
                    rng
                )
        elseif simtype == :fixedtime || simtype == "fixedtime"
            population = 
                simulatefixedtime!(
                    population, 
                    input.tmax, 
                    input.maxmodules, 
                    input.b, 
                    input.d, 
                    input.bdrate, 
                    input.branchrate, 
                    input.modulesize, 
                    input.branchinitsize, 
                    rng
                )
        else
            error("incorrect value for simtype")
        end
        if length(population) != 0
            break
        end
    end
    population = 
        processresults!(population, input.μ, input.clonalmutations, rng)
    return MultiSimulation(input, population)
end

function runsimulation(::Type{Cell}, input::MultilevelBranchingMoranInput, rng::AbstractRNG=Random.GLOBAL_RNG) 
    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            input.modulesize, 
            clonalmutations=0
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
                1, #assume a single mutation at each division and add distribution of
                :fixed, #mutations later
                rng,
                moduleupdate=:moran
            )
        if length(population) != 0
            break
        end
    end
    population = 
        processresults!(population, input.μ, input.clonalmutations, rng)
    return MultiSimulation(input, population)
end

function runsimulation_timeseries_returnfinalpop(::Type{Cell}, input::MultilevelBranchingMoranInput, timesteps, func, rng::AbstractRNG=Random.GLOBAL_RNG) 
    population = initialize_population(
        input.modulesize, 
        clonalmutations=0
    )
    data = []
    t0 = 0.0
    for t in timesteps
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
            rng;
            moduleupdate=:moran,
            t0
        )
        length(population) < input.maxmodules || break
        push!(data, func(population))
        t0 = t
    end
    return data, population
end

function runsimulation_timeseries(::Type{T}, input::MultilevelBranchingMoranInput, timesteps, func, rng::AbstractRNG=Random.GLOBAL_RNG) where T
    return runsimulation_timeseries_returnfinalpop(T, input, timesteps, func, rng)[1] 
end
"""
    simulate!(population, tmax, maxmodules, b, d, bdrate, branchrate, 
        modulesize, branchinitsize, rng; moduleupdate=:branching[, t0])

Run a single multilevel simulation using the Gillespie algorithm, with the 
`population` giving the initial state. Simulation runs until the module population 
size reaches `maxmodules` or the age of the population reaches `tmax`. 

"""
function simulate!(population, tmax, maxmodules, b, d, bdrate, branchrate, 
    modulesize, branchinitsize, μ, mutationdist, rng; moduleupdate=:branching, t0=nothing)

    nextID = getnextID(population) #get next id
    nextmoduleID = maximum(cellmodule.id for cellmodule in population) + 1
    t = isnothing(t0) ? age(population) : t0
    while t < tmax && (moduleupdate==:moran || length(population) < maxmodules)
        population, t, nextID = 
            update_population!(population, b, d, bdrate, branchrate, modulesize, 
                branchinitsize, t, nextID, nextmoduleID, μ, mutationdist, tmax, maxmodules, rng; moduleupdate)
        #returns empty list of modules if population dies out
        if length(population) == 0
            return population
        end
    end
    return population
end

"""
    update_population!(population, b, d, bdrate, branchrate, modulesize, 
        branchinitsize, nextID, t, tmax, rng)

Perform a single Gillespie step to advance the simulation of a multilevel population, and
increase the time t → t + Δt.

Calculates the transition rates for a Moran, birth, death or module branching event; then 
selects from these transitions with probability proportional to rate. Time is increased by
`Δt ~ Exp(sum(transitionrates))`. If new time exceeds `tmax` do not perform any transition.
"""

function update_population!(population, b, d, bdrate, branchrate, modulesize, 
    branchinitsize, t, nextID, nextmoduleID, μ, mutationdist, tmax, maxmodules, rng; moduleupdate=:branching)

    transitionrates = get_transitionrates(population, b, d, 
        bdrate, branchrate, modulesize)
    t += exptime(rng, sum(transitionrates))
    #only update the population if t < tmax
    if t < tmax
        #choose transition type: 1=moran, 2=birth, 3=death, 4=branch
        transitionid = sample(
            rng, 
            1:4, 
            ProbabilityWeights(transitionrates ./ sum(transitionrates))
        )
        population, nextID, nextmoduleID = transition!(
            population, 
            transitionid, 
            modulesize, 
            branchinitsize, 
            t, 
            nextID,
            nextmoduleID,
            μ, 
            mutationdist,
            maxmodules,
            rng;
            moduleupdate
        )
    end
    return population, t, nextID, nextmoduleID
end

"""
    transition!(population, transitionid, modulesize, branchinitsize, nextID, t, 
        μ, rng)
    
Perform a single transition step on `population`, determined by `transitionid`.
"""
function transition!(population, transitionid, modulesize, branchinitsize, t, nextID, nextmoduleID,
    μ, mutationdist, maxmodules, rng; moduleupdate=:branching)
    
    if transitionid == 1
        _, nextID = moranupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 2
        _, nextID = birthupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 3
        deathupdate!(population, modulesize, t, μ, mutationdist, rng)
    elseif transitionid == 4
        if moduleupdate == :branching || length(population) < maxmodules
            _, nextmoduleID = modulebranchingupdate!(population, nextmoduleID, modulesize, branchinitsize, t, rng)
        elseif moduleupdate == :moran
            _, nextmoduleID = modulemoranupdate!(population, nextmoduleID, modulesize, branchinitsize, t, μ, mutationdist, rng)
        end
    end
    return population, nextID, nextmoduleID
end

"""
    moranupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a homeostatic module uniformly at random to undergo a single Moran update. From 
that module one cell divides and one cell dies.
"""
function moranupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    cellmodule, parentcell, deadcell = 
        choose_homeostaticmodule_cells(population, modulesize, rng)
    _, nextID = celldivision!(cellmodule, parentcell, t, nextID, μ, mutationdist, rng)
    celldeath!(cellmodule, deadcell, t, μ, mutationdist, rng)
    updatemodulehistory!(cellmodule, 0, t)
    return population, nextID
end

"""
    birthupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    cellmodule, parentcell = 
        choose_growingmodule_cell(population, modulesize, rng)
    _, nextID = celldivision!(cellmodule, parentcell, t, nextID, μ, mutationdist, rng)
    updatemodulehistory!(cellmodule, 1, t)
    return population, nextID
end


"""
    deathupdate!(population, modulesize, t, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to die. If 
cell death results in an empty module, remove that module from the population.
"""

function deathupdate!(population, modulesize, t, μ, mutationdist, rng)
    cellmodule, deadcell = 
        choose_growingmodule_cell(population, modulesize, rng)
    celldeath!(cellmodule, deadcell, t, μ, mutationdist, rng)
    updatemodulehistory!(cellmodule, -1, t)
    if length(cellmodule) == 0
        moduledeath!(population, cellmodule, t, μ, mutationdist, rng)
    end
    return population
end

"""
    modulebranchingupdate!(population, modulesize, branchinitsize, t, rng)
Select a homeostatic module, uniformly at random, to undergo branching. Cells are sampled 
(number of cells giving by `branchinitsize`) from the parent module to form a new module, 
which is added to `population`
"""
function modulebranchingupdate!(population, nextmoduleID, modulesize, branchinitsize, t, rng)
    parentmodule = choose_homeostaticmodule(population, modulesize, rng)
    _, newmodule, nextmoduleID = 
        sample_new_module!(parentmodule, nextmoduleID, branchinitsize, t, rng)
    push!(population, newmodule)
    return population, nextmoduleID
end

function modulemoranupdate!(population, nextmoduleID, modulesize, branchinitsize, t, μ, mutationdist, rng)
    parentmodule = choose_homeostaticmodule(population, modulesize, rng)
    _, newmodule, nextmoduleID = 
        sample_new_module!(parentmodule, nextmoduleID, branchinitsize, t, rng)
    push!(population, newmodule)
    deadmodule = rand(rng, population)
    moduledeath!(population, deadmodule, t, μ, mutationdist, rng)
    return population, nextmoduleID
end

"""
    updatemodulehistory!(cellmodule, ΔN, newtime)
Adds new population size (old size + `ΔN`) and `newtime` to the end of `cellmodule.Nvec`
and `cellmodule.tvec`.

"""
function updatemodulehistory!(cellmodule, ΔN, newtime)
    push!(cellmodule.Nvec, cellmodule.Nvec[end] + ΔN)
    push!(cellmodule.tvec, newtime)
    return cellmodule
end

"""
    choose_homeostaticmodule(population, maxmodulesize, rng::AbstractRNG)
Select a homeostatic module, i.e. module of size `maxmodulesize`, uniformly at random.
"""
function choose_homeostaticmodule(population, maxmodulesize, rng::AbstractRNG)
    homeostatic_modules = filter(x -> length(x) == maxmodulesize, population)
    return rand(rng, homeostatic_modules)
end

"""
    choose_homeostaticmodule_cells(population, maxmodulesize, rng::AbstractRNG)
Select a homeostatic module and the ids of two cells from the module, uniformly at random.
"""
function choose_homeostaticmodule_cells(population, maxmodulesize, rng::AbstractRNG)
    chosenmodule = choose_homeostaticmodule(population, maxmodulesize, rng)
    return chosenmodule, rand(rng, 1:maxmodulesize), rand(rng, 1:maxmodulesize)
end

"""
    choose_growingmodule_cells(population, maxmodulesize, rng::AbstractRNG)
Select a cell uniformly at random from all cells in growing (non-homeostatic) modules, and
return the module and cell id.
"""
function choose_growingmodule_cell(population, maxmodulesize, rng::AbstractRNG)
    modulesizes = map(length, population)
    modules = population[modulesizes .< maxmodulesize]
    modulesizes = modulesizes[modulesizes .< maxmodulesize]
    chosenmodule = sample(rng, modules, ProbabilityWeights(modulesizes ./ sum(modulesizes)))
    chosencell = rand(rng, 1:length(chosenmodule))
    return chosenmodule, chosencell
end

"""
    get_transitionrates(population, b, d, bdrate, branchrate, modulesize)
Compute the rates for moran update, cell division, death and module branching and return as 
a Vector.
"""
function get_transitionrates(population, b, d, bdrate, branchrate, modulesize)
    rates = zeros(Float64, 4)
    number_homeostatic_modules = 0
    for (i, cellmodule) in enumerate(population)
        N = length(cellmodule)
        if N < modulesize
            rates[2] += N * b
            rates[3] += N * d
        elseif N == modulesize
            number_homeostatic_modules += 1
        else 
            error("module size exceeds homeostatic size")
        end
    end
    rates[1] = number_homeostatic_modules * bdrate * modulesize
    rates[4] = number_homeostatic_modules * branchrate
    return rates
end

"""
    moduledeath!(population, cellmodule)
Kill all live cells in `cellmodule` and remove it from `population`.
    """
function moduledeath!(population, cellmodule, t, μ, mutationdist, rng)
    length(cellmodule) == 0 || killallcells!(cellmodule.cells, t, μ, mutationdist, rng)
    return deleteat!(population, findall(x->x==cellmodule, population))
end

function moduledeath!(population, cellmodule)
    @assert length(cellmodule) == 0 "trying to delete non-empty module"
    return deleteat!(population, findall(x->x==cellmodule, population))
end

killallcells!(population::Vector{Cell}, args...) = nothing

"""
    simulatefixedtime!(population, tmax, maxmodules, b, d, bdrate, branchrate, 
        modulesize, branchinitsize, rng)

Run a single multilevel simulation until `tmax`, with the `population` giving the 
initial state. Faster than `simulate!`, but cannot run to fixed population size.
"""
function simulatefixedtime!(population, tmax, maxmodules, b, d, bdrate, branchrate, 
    modulesize, branchinitsize, rng)
    
    for cellmodule in population
        check_module_number(population, maxmodules)
        while true
            cellmodule, newcellmodule =
                module_simulate_to_branching!(
                    cellmodule, 
                    tmax, 
                    b, 
                    d, 
                    bdrate, 
                    branchrate, 
                    modulesize, 
                    branchinitsize,
                    length(population) + 1,
                    rng
                )
                
            if newcellmodule !== nothing
                push!(population, newcellmodule)
            elseif length(cellmodule) == 0
                moduledeath!(population, cellmodule)
                break
            else
                break
            end
        end
    end
    return population
end

"""
    check_module_number(population, maxmodules)
Throws an error if the number of modules in `population` exceeds `maxmodules`.
"""
function check_module_number(population, maxmodules)
    numbermodules = length(population)    
    if numbermodules >= maxmodules
        error("population size exceeds maxmodules: ", numbermodules, " > ", maxmodules)
    end
end

"""
    module_simulate_to_branching(cellmodule, tmax, b, d, bdrate, branchrate, 
        modulesize, branchinitsize, newmoduleid, rng)

    Simulates cellmodule dynamics until a module branching event occurs, or tmax is
    reached.

    If cellmodule is smaller than `modulesize` it grows according to a branching 
    process, with birthrate `b` and deathrate `d`, until it reaaches `modukesize`. Once it 
    reaches `modulesize` the dynamic switches to a Moran process with event rate `bdrate`. 

"""

function module_simulate_to_branching!(cellmodule, tmax, b, d, bdrate, 
    branchrate, modulesize, branchinitsize, newmoduleid, rng)

    #branching process until module reaches homeostatic size
    if length(cellmodule) < modulesize
        branchingprocess!(
            cellmodule, 
            b, 
            d, 
            modulesize, 
            1, 
            :fixed,
            tmax,
            rng; 
            numclones=0, 
            maxclonesize=Inf
        )

        if length(cellmodule) == 0
            return cellmodule, nothing
        end
    end

    #get time of next branching event
    branchtime = 1 / branchrate .* exptime(rng) + cellmodule.tvec[end]

    #if module branching will occur before maxime is reached, simulate moran process until
    #branchtime and sample new module
    if branchtime < tmax
        moranprocess!(cellmodule, bdrate, branchtime, 1, :fixed, rng)
        cellmodule, newcellmodule = 
            sample_new_module!(cellmodule, newmoduleid, branchinitsize, branchtime, rng)
    
        return cellmodule, newcellmodule

    #if module branching will occur after tmax is reached, simulate moran process until
    #tmax 
    else
        moranprocess!(cellmodule, bdrate, tmax, 1, :fixed, rng)
        return cellmodule, nothing
    end
end

"""
    sample_new_module(cellmodule, newmoduleid, branchinitsize, branchtime, 
        rng::AbstractRNG)

Sample cells, uniformly at random, from `cellmodule` to form `newcellmodule`. The
number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module!(cellmodule::T, nextmoduleID, branchinitsize, branchtime, 
    rng::AbstractRNG) where T<: AbstractModule

    sampleids = sample(rng, 1:length(cellmodule.cells), branchinitsize, replace=false)
    newcellmodule = 
        initialize_from_cells(T, cellmodule.cells[sampleids], cellmodule.subclones, 
            nextmoduleID, cellmodule.id, inittime=branchtime)
    cellremoval!(cellmodule, sampleids)
    push!(cellmodule.Nvec, cellmodule.Nvec[end] - branchinitsize)
    push!(cellmodule.tvec, branchtime)

    return cellmodule, newcellmodule, nextmoduleID
end

"""
    initialize_population([modulesize=nothing]; clonalmutations=0)

Create a Vector{CellModule} of length 1 that contains a single module comprising a single 
cell at time t=0. 
"""
function initialize_population(modulesize=nothing; clonalmutations=0)

    #population is defined by CellModule vector with each entry coprresponding to
    #an individual module
    population = CellModule[]

    #initialise population with a single module conisting of a single cell
    initialmodule = initializesim_branching(
        modulesize,
        clonalmutations=clonalmutations
    )
    push!(population, initialmodule)

    return population
end

getnextID(population::Vector{CellModule}) = maximum(map(x->getnextID(x.cells), population))