function run_multilevel_from_file(filename, outputdir)
    inputdict = JSON.parsefile(filename, dicttype=Dict{Symbol, Any})
    input = loadinput(MultilevelInput, inputdict[:input])
    seeds = inputdict[:seeds]
    output = inputdict[:output]
    if inputdict[:repeat] > 1
        for i in 1:inputdict[:repeat]
            rng, seed = seeds !== nothing ? (MersenneTwister(seeds[i]), seeds[i]) : (MersenneTwister(), nothing)
            population = multilevel_simulation(input, rng, Symbol(inputdict[:simtype]))
            saveoutput(population, output, outputdir, seed, i)
        end
    else
        rng = MersenneTwister(seeds)
        population = multilevel_simulation(input, rng, Symbol(inputdict[:simtype]))
        saveoutput(population, output, outputdir, seeds, 1)
    end
end

function run_multilevel_from_file(filename, outputdir, id)
    inputdict = JSON.parsefile(filename, dicttype=Dict{Symbol, Any})
    input = loadinput(MultilevelInput, inputdict[:input])
    seed = (inputdict[:seeds])[id]
    output = inputdict[:output]
    rng = MersenneTwister(seed)
    population = multilevel_simulation(input, rng, Symbol(inputdict[:simtype]))
    saveoutput(population, output, outputdir, seed, id)
end

"""
    multilevel_simulation(input::MultilevelInput, rng::AbstractRNG = MersenneTwister; 
        <keyword arguments>)

    Simulate a growing population of modules, until it reaches `maxmodules`.
    
    Start with a single module comprised of a single cell that grows via a branching process  
    to size `input.modulesize`, then switches to a Moran process. 
    New modules are created with rate `input.branchrate`, by sampling cells from the parent 
    module. Return output as a Population.

    - if `simtype == :normal`, simulation is
    implemented by a Gillespie algorithm and runs until the number of modules exceeds 
    input.maxmodules or time exceeds input.maxtime.

    - if `simtype == :fixedtime`, a list of modules is created (initially of length 1)
    and each module is simulated independently until input.maxtime is reached. New modules 
    are added to the end of the list. This implementation is faster, but only works if we 
    end the simulation at a fixed time (not fixed size). If the list of modules is bigger 
    than input.maxmodules and error is thrown.

"""
function multilevel_simulation(input::MultilevelInput, rng::AbstractRNG=Random.GLOBAL_RNG, 
    simtype=:normal)

    #if the population dies out we start a new simulation
    while true 
        populationtracker = initialize_population(
            input.modulesize, 
            clonalmutations=0
        )
        if simtype == :normal || simtype == "normal"
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
                    rng
                )
        elseif simtype == :fixedtime || simtype == "fixedtime"
            populationtracker = 
                simulatefixedtime!(
                    populationtracker, 
                    input.maxtime, 
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
        if length(populationtracker) != 0
            break
        end
    end
    populationtracker = 
        processresults!(populationtracker, input.μ, input.clonalmutations, rng)
    return Population(input, populationtracker)
end

"""
    simulate!(populationtracker, maxtime, maxmodules, b, d, bdrate, branchrate, 
        modulesize, branchinitsize, rng)

Run a single multilevel simulation using the Gillespie algorithm, with the 
`populationtracker` giving the initial state. Simulation runs until the module population 
size reaches `maxmodules` or the age of the population reaches `maxtime`. 

"""
function simulate!(populationtracker, maxtime, maxmodules, b, d, bdrate, branchrate, 
    modulesize, branchinitsize, rng)

    mutID = getmutID(populationtracker[1].cells) #get id of next mutation
    t = age(populationtracker)
    while t < maxtime && length(populationtracker) < maxmodules
        populationtracker, mutID, t = 
            update_population!(populationtracker, b, d, bdrate, branchrate, modulesize, 
                branchinitsize, mutID, t, maxtime, rng)
        #returns empty list of modules if population dies out
        if length(populationtracker) == 0
            return populationtracker
        end
    end
    return populationtracker
end

"""
    update_population!(populationtracker, b, d, bdrate, branchrate, modulesize, 
        branchinitsize, mutID, t, maxtime, rng)

Perform a single Gillespie step to advance the simulation of a multilevel population, and
increase the time t → t + Δt.

Calculates the transition rates for a Moran, birth, death or module branching event; then 
selects from these transitions with probability proportional to rate. Time is increased by
`Δt ~ Exp(sum(transitionrates))`. If new time exceeds `maxtime` do not perform any transition.
"""

function update_population!(populationtracker, b, d, bdrate, branchrate, modulesize, 
    branchinitsize, mutID, t, maxtime, rng)

    transitionrates = get_transitionrates(populationtracker, b, d, 
        bdrate, branchrate, modulesize)
    t += exptime(rng, sum(transitionrates))
    #only update the population if t < maxtime
    if t < maxtime
        #choose transition type: 1=moran, 2=birth, 3=death, 4=branch
        transitionid = sample(
            rng, 
            1:4, 
            ProbabilityWeights(transitionrates ./ sum(transitionrates))
        )
        populationtracker, mutID = transition!(
            populationtracker, 
            transitionid, 
            modulesize, 
            branchinitsize, 
            mutID, 
            t, 
            1, 
            rng
        )
    end
    return populationtracker, mutID, t
end

"""
    transition!(populationtracker, transitionid, modulesize, branchinitsize, mutID, t, 
        μ, rng)
    
Perform a single transition step on `populationtracker`, determined by `transitionid`.
"""
function transition!(populationtracker, transitionid, modulesize, branchinitsize, mutID, t, 
    μ, rng)
    
    if transitionid == 1
        _, mutID = moranupdate!(populationtracker, modulesize, mutID, t, μ, rng)
    elseif transitionid == 2
        _, mutID = birthupdate!(populationtracker, modulesize, mutID, t, μ, rng)
    elseif transitionid == 3
        deathupdate!(populationtracker, modulesize, t, rng)
    elseif transitionid == 4
        modulebranchingupdate!(populationtracker, modulesize, branchinitsize, t, rng)
    end
    return populationtracker, mutID
end

"""
    moranupdate!(populationtracker, modulesize, mutID, t, μ, rng)

Selects a homeostatic module uniformly at random to undergo a single Moran update. From 
that module one cell divides and one cell dies.
"""
function moranupdate!(populationtracker, modulesize, mutID, t, μ, rng)
    moduletracker, parentcell, deadcell = 
        choose_homeostaticmodule_cells(populationtracker, modulesize, rng)
    _, mutID = celldivision!(moduletracker, parentcell, mutID, μ, rng)
    celldeath!(moduletracker, deadcell)
    updatemodulehistory!(moduletracker, 0, t)
    return populationtracker, mutID
end

"""
    birthupdate!(populationtracker, modulesize, mutID, t, μ, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate!(populationtracker, modulesize, mutID, t, μ, rng)
    moduletracker, parentcell = 
        choose_growingmodule_cell(populationtracker, modulesize, rng)
    _, mutID = celldivision!(moduletracker, parentcell, mutID, μ, rng, fixedmu=true)
    updatemodulehistory!(moduletracker, 1, t)
    return populationtracker, mutID
end

"""
    deathupdate!(populationtracker, modulesize, t, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to die. If 
cell death results in an empty module, remove that module from the population.
"""
function deathupdate!(populationtracker, modulesize, t, rng)
    moduletracker, deadcell = 
        choose_growingmodule_cell(populationtracker, modulesize, rng)
    celldeath!(moduletracker, deadcell)
    updatemodulehistory!(moduletracker, -1, t)
    if length(moduletracker) == 0
        moduledeath!(populationtracker, moduletracker)
    end
    return populationtracker
end

"""
    modulebranchingupdate!(populationtracker, modulesize, branchinitsize, t, rng)
Select a homeostatic module, uniformly at random, to undergo branching. Cells are sampled 
(number of cells giving by `branchinitsize`) from the parent module to form a new module, 
which is added to `populationtracker`
"""
function modulebranchingupdate!(populationtracker, modulesize, branchinitsize, t, rng)
    parentmodule = choose_homeostaticmodule(populationtracker, modulesize, rng)
    _, newmodule = 
        sample_new_module!(parentmodule, length(populationtracker) + 1, 
            branchinitsize, t, rng)
    push!(populationtracker, newmodule)
    return populationtracker
end

"""
    updatemodulehistory!(moduletracker, ΔN, newtime)
Adds new population size (old size + `ΔN`) and `newtime` to the end of `moduletracker.Nvec`
and `moduletracker.tvec`.

"""
function updatemodulehistory!(moduletracker, ΔN, newtime)
    push!(moduletracker.Nvec, moduletracker.Nvec[end] + ΔN)
    push!(moduletracker.tvec, newtime)
    return moduletracker
end

"""
    choose_homeostaticmodule(populationtracker, maxmodulesize, rng::AbstractRNG)
Select a homeostatic module, i.e. module of size `maxmodulesize`, uniformly at random.
"""
function choose_homeostaticmodule(populationtracker, maxmodulesize, rng::AbstractRNG)
    homeostatic_modules = filter(x -> length(x) == maxmodulesize, populationtracker)
    return rand(rng, homeostatic_modules)
end

"""
    choose_homeostaticmodule_cells(populationtracker, maxmodulesize, rng::AbstractRNG)
Select a homeostatic module and two cells from the module, uniformly at random.
"""
function choose_homeostaticmodule_cells(populationtracker, maxmodulesize, rng::AbstractRNG)
    chosenmodule = choose_homeostaticmodule(populationtracker, maxmodulesize, rng)
    return chosenmodule, rand(rng, 1:maxmodulesize), rand(rng, 1:maxmodulesize)
end

"""
    choose_growingmodule_cells(populationtracker, maxmodulesize, rng::AbstractRNG)
Select a cell unifor´mly at random from all cells in growing (non-homeostatic) modules, and
return the module and cell.
"""
function choose_growingmodule_cell(populationtracker, maxmodulesize, rng::AbstractRNG)
    modulesizes = map(length, populationtracker)
    modules = populationtracker[modulesizes .< maxmodulesize]
    modulesizes = modulesizes[modulesizes .< maxmodulesize]
    chosenmodule = sample(rng, modules, ProbabilityWeights(modulesizes ./ sum(modulesizes)))
    chosencell = rand(rng, 1:length(chosenmodule))
    return chosenmodule, chosencell
end

"""
    get_transitionrates(populationtracker, b, d, bdrate, branchrate, modulesize)
Compute the rates for moran update, cell division, death and module branching and return as 
a Vector.
"""
function get_transitionrates(populationtracker, b, d, bdrate, branchrate, modulesize)
    rates = zeros(Float64, 4)
    number_homeostatic_modules = 0
    for (i, moduletracker) in enumerate(populationtracker)
        N = length(moduletracker)
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
    moduledeath!(populationtracker, moduletracker)
Remove the given `moduletracker` from `populationtracker`.
"""
function moduledeath!(populationtracker, moduletracker)
    return deleteat!(populationtracker, findall(x->x==moduletracker, populationtracker))
end

"""
    simulatefixedtime!(populationtracker, maxtime, maxmodules, b, d, bdrate, branchrate, 
        modulesize, branchinitsize, rng)

Run a single multilevel simulation until `maxtime`, with the `populationtracker` giving the 
initial state. Faster than `simulate!`, but cannot run to fixed population size.
"""
function simulatefixedtime!(populationtracker, maxtime, maxmodules, b, d, bdrate, branchrate, 
    modulesize, branchinitsize, rng)
    
    for moduletracker in populationtracker
        check_module_number(populationtracker, maxmodules)
        while true
            moduletracker, newmoduletracker =
                module_simulate_to_branching!(
                    moduletracker, 
                    maxtime, 
                    b, 
                    d, 
                    bdrate, 
                    branchrate, 
                    modulesize, 
                    branchinitsize,
                    length(populationtracker) + 1,
                    rng
                )
                
            if newmoduletracker !== nothing
                push!(populationtracker, newmoduletracker)
            elseif length(moduletracker) == 0
                moduledeath!(populationtracker, moduletracker)
                break
            else
                break
            end
        end
    end
    return populationtracker
end

"""
    check_module_number(populationtracker, maxmodules)
Throws an error if the number of modules in `populationtracker` exceeds `maxmodules`.
"""
function check_module_number(populationtracker, maxmodules)
    numbermodules = length(populationtracker)    
    if numbermodules >= maxmodules
        error("population size exceeds maxmodules: ", numbermodules, " > ", maxmodules)
    end
end

"""
    module_simulate_to_branching(moduletracker, maxtime, b, d, bdrate, branchrate, 
        modulesize, branchinitsize, newmoduleid, rng)

    Simulates moduletracker dynamics until a module branching event occurs, or maxtime is
    reached.

    If moduletracker is smaller than `modulesize` it grows according to a branching 
    process, with birthrate `b` and deathrate `d`, until it reaaches `modukesize`. Once it 
    reaches `modulesize` the dynamic switches to a Moran process with event rate `bdrate`. 

"""

function module_simulate_to_branching!(moduletracker, maxtime, b, d, bdrate, 
    branchrate, modulesize, branchinitsize, newmoduleid, rng)

    #branching process until module reaches homeostatic size
    if length(moduletracker) < modulesize
        moduletracker = 
            branchingprocess!(moduletracker, b, d, modulesize, 1, 
                rng, numclones=0, fixedmu=true, maxclonesize=Inf, tmax=maxtime)
        if length(moduletracker) == 0
            return moduletracker, nothing
        end
    end

    #get time of next branching event
    branchtime = 1 / branchrate .* exptime(rng) + moduletracker.tvec[end]

    #if module branching will occur before maxime is reached, simulate moran process until
    #branchtime and sample new module
    if branchtime < maxtime
        moduletracker = 
            moranprocess!(moduletracker, bdrate, branchtime, 1, rng, numclones=0, 
                fixedmu=true)
        
        moduletracker, newmoduletracker = 
            sample_new_module!(moduletracker, newmoduleid, branchinitsize, branchtime, rng)
    
        return moduletracker, newmoduletracker

    #if module branching will occur after maxime is reached, simulate moran process until
    #maxtime 
    else
        moduletracker = 
            moranprocess!(moduletracker, bdrate, maxtime, 1, rng, numclones=0, 
                fixedmu=true)
        return moduletracker, nothing
    end
end

"""
    sample_new_module(moduletracker, newmoduleid, branchinitsize, branchtime, 
        rng::AbstractRNG)

Sample cells, uniformly at random, from `moduletracker` to form `newmoduletracker`. The
number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module!(moduletracker, newmoduleid, branchinitsize, branchtime, 
    rng::AbstractRNG)

    sampleids = sample(rng, 1:length(moduletracker.cells), branchinitsize, replace=false)
    newmoduletracker = 
        initializesim_from_cells(moduletracker.cells[sampleids], moduletracker.subclones, 
            newmoduleid, moduletracker.id, inittime=branchtime)
    moduletracker = celldeath!(moduletracker, sampleids)
    push!(moduletracker.Nvec, moduletracker.Nvec[end] - branchinitsize)
    push!(moduletracker.tvec, branchtime)

    return moduletracker, newmoduletracker
end

"""
    initialize_population([modulesize=nothing]; clonalmutations=0)

Create a Vector{ModuleTracker} of length 1 that contains a single module comprising a single 
cell at time t=0. 
"""
function initialize_population(modulesize=nothing; clonalmutations=0)

    #population is defined by ModuleTracker vector with each entry coprresponding to
    #an individual module
    population = ModuleTracker[]

    #initialise population with a single module conisting of a single cell
    initialmodule = initializesim_branching(
        modulesize,
        clonalmutations=clonalmutations
    )
    push!(population, initialmodule)

    return population
end