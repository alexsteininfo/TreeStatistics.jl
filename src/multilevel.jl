"""
    multilevel_simulation(input::MultilevelInput, rng::AbstractRNG = MersenneTwister; 
        <keyword arguments>)

    Run a simulation process that models somatic evolution in a growing population. Starts 
    with a single cell which grows via a branching process to a "module" of specified size, 
    then switches to a Moran process. New modules are created with given rate, by 
    sampling cells from the parent module.
"""

function multilevel_simulation(input::MultilevelInput; 
    rng::AbstractRNG=Random.GLOBAL_RNG, maxmodules=1e6, showprogress=false)

    populationtracker = initialize_population(
        input.modulesize, 
        clonalmutations=input.clonalmutations
        )

    mutID = getmutID(populationtracker[1].cells)
    t = 0

    while t < input.pop_age && length(populationtracker) < maxmodules
    
        transitionrates = get_transitionrates(populationtracker, input.b, input.d, 
            input.bdrate, input.branchrate, input.modulesize)
        t += 1/sum(transitionrates) .* exptime(rng)
        # r = rand(rng)
        # 1=moran, 2=birth, 3=death, 4=branch
        transitionid = sample(rng, 1:4, ProbabilityWeights(transitionrates ./ sum(transitionrates)))
        if transitionid == 1
            #moran update of a homeostatic module
            moduletracker, parentcell, deadcell = 
                choose_homeostaticmodule_cells(populationtracker, input.modulesize, rng)
            _, mutID = celldivision!(moduletracker, parentcell, mutID, 1, rng)
            celldeath!(moduletracker, deadcell)
            update!(moduletracker, 0, t)

        elseif transitionid == 2
            #cell division in a growing module
            moduletracker, parentcell = 
                choose_growingmodule_cell(populationtracker, input.modulesize, rng)
            _, mutID = celldivision!(moduletracker, parentcell, mutID, 1, rng, fixedmu=true)
            update!(moduletracker, 1, t)

        elseif transitionid == 3
            #cell division in a growing module
            _, deadcell = 
                choose_nonhomeostatic_cell(populationtracker, input.modulesize, rng)
            celldeath!(moduletracker, deadcell)
            update!(moduletracker, -1, t)

        elseif transitionid == 4
            #branching of a homeostatic module to create a new module
            parentmodule = choose_homeostaticmodule(populationtracker, input.modulesize, rng)
            _, newmodule = 
                sample_new_module!(parentmodule, length(populationtracker) + 1, 
                    input.branchinitsize, t, rng)
            push!(populationtracker, newmodule)
        end
    end
    populationtracker = processresults!(populationtracker, input.μ, input.clonalmutations, rng)
    return Population(input, populationtracker)
end

function update!(moduletracker, ΔN, newtime)
    push!(moduletracker.Nvec, moduletracker.Nvec[end] + ΔN)
    push!(moduletracker.tvec, newtime)
    return moduletracker
end


function choose_homeostaticmodule(populationtracker, maxmodulesize, rng::AbstractRNG)
    homeostatic_modules = filter(x -> length(x) == maxmodulesize, populationtracker)
    return rand(rng, homeostatic_modules)
end

function choose_homeostaticmodule_cells(populationtracker, maxmodulesize, rng::AbstractRNG)
    chosenmodule = choose_homeostaticmodule(populationtracker, maxmodulesize, rng)
    return chosenmodule, rand(rng, 1:maxmodulesize), rand(rng, 1:maxmodulesize)
end

function choose_growingmodule_cell(populationtracker, maxmodulesize, rng::AbstractRNG)
    modulesizes = map(length, populationtracker)
    modules = populationtracker[modulesizes .< maxmodulesize]
    modulesizes = modulesizes[modulesizes .< maxmodulesize]
    chosenmodule = sample(rng, modules, ProbabilityWeights(modulesizes ./ sum(modulesizes)))
    chosencell = rand(rng, 1:length(chosenmodule))
    return chosenmodule, chosencell
end

function get_transitionrates(populationtracker, b, d, bdrate, branchrate, modulesize)
    #vector of rates for moran update, cell division, death and module branching
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

function multilevel_simulation_fast(input::MultilevelInput; 
    rng::AbstractRNG=Random.GLOBAL_RNG, maxmodules=1e6, showprogress=false)

    populationtracker = initialize_population(
        input.modulesize, 
        clonalmutations=input.clonalmutations
        )
    
    for moduletracker in populationtracker
        check_module_number(moduletracker.id, populationtracker, maxmodules, showprogress)
        while true
            moduletracker, newmoduletracker =
                module_simulate_to_branching!(
                    moduletracker, 
                    input,  
                    length(populationtracker) + 1,
                    rng
                )
            if newmoduletracker !== nothing
                push!(populationtracker, newmoduletracker)
            else
                break
            end
        end
    end
    populationtracker = processresults!(populationtracker, input.μ, input.clonalmutations, rng)
    return Population(input, populationtracker)
end

function check_module_number(moduleid, populationtracker, maxmodules, showprogress)
    numbermodules = length(populationtracker)    
    if numbermodules >= maxmodules
        error("population size exceeds maxmodules: ", numbermodules, " > ", maxmodules)
    end
    if showprogress
        print("\rmodule ",moduleid, " of ", numbermodules, "\t")
    end
end

function module_simulate_to_branching!(moduletracker::ModuleTracker, input::MultilevelInput, 
    newmoduleid, rng::AbstractRNG=Random.GLOBAL_RNG
    )

    if moduletracker.Nvec[end] < input.modulesize
        moduletracker = 
            branchingprocess!(moduletracker, input.b, input.d, input.modulesize, 1, 
                rng, numclones=input.numclones, fixedmu=true, selection=input.selection, 
                tevent=input.tevent, maxclonesize=Inf, tmax=input.pop_age)
    end 
    branchtime = 1 / input.branchrate .* exptime(rng) + moduletracker.tvec[end]
    if branchtime < input.pop_age
        moduletracker = 
            moranprocess!(moduletracker, input.bdrate, branchtime, 1, rng::AbstractRNG; 
            numclones = input.numclones, fixedmu = true, selection = input.selection, 
            tevent = input.tevent)
        
        moduletracker, newmoduletracker = 
            sample_new_module!(moduletracker, newmoduleid, input.branchinitsize, branchtime, rng)
    
        return moduletracker, newmoduletracker
    else
        moduletracker = 
            moranprocess!(moduletracker, input.bdrate, input.pop_age, 1, rng::AbstractRNG; 
            numclones = input.numclones, fixedmu = true, selection = input.selection, 
            tevent = input.tevent)
        return moduletracker, nothing
    end
end

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
