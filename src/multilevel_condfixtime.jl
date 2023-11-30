function runsimulation_condfixtime(input::MultilevelInput, args...)
    return runsimulation_condfixtime(Cell, WellMixed, input, args...) 
end

function runsimulation_condfixtime(::Type{Cell}, ::Type{S}, input::MultilevelBranchingInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG, simtype=:normal) where S <: ModuleStructure
    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            T, S, 
            clonalmutations, 
            getNinit(input),
            getmaxmodulesize(input),
            input.birthrate,
            input.deathrate,
            input.moranrate,
            input.asymmetricrate;
            rng
        )
        nextID, nextmoduleID = 2, 2

        population, condfixtimes = 
            simulate_condfixtime!(
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
                rng
            )
        if length(population) != 0
            return condfixtimes, MultiSimulation(input, population)
        end
    end
end

function runsimulation_condfixtime(::Type{Cell}, ::Type{S}, 
    input::MultilevelBranchingMoranInput, rng::AbstractRNG=Random.GLOBAL_RNG)  where S <: ModuleStructure
    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            Cell,
            S,
            clonalmutations=input.clonalmutations
        )

        nextID, nextmoduleID = 2, 2
        population, condfixtimes = 
            simulate_condfixtime!(
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
                moduleupdate=:moran
            )
        if length(population) != 0
            return condfixtimes, MultiSimulation(input, population)
        end
    end
end

function runsimulation_condfixtime_to_nfixed(::Type{Cell}, ::Type{S}, 
    input::MultilevelBranchingMoranInput, nfixed, rng::AbstractRNG=Random.GLOBAL_RNG) where S <: ModuleStructure
    #if the population dies out we start a new simulation
    while true 
        population = initialize_population(
            Cell,
            S,
            input.modulesize, 
            clonalmutations=input.clonalmutations
        )

        nextID, nextmoduleID = 2, 2
        population, condfixtimes = 
            simulate_condfixtime_to_nfixed!(
                population, 
                nfixed,
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
                moduleupdate=:moran
            )
        if length(population) != 0
            return condfixtimes, MultiSimulation(input, population)
        end
    end
end


"""
    simulate_condfixtime!(population, tmax, maxmodules, birthrate, deathrate, moranrate, branchrate, 
        modulesize, branchinitsize, rng; moduleupdate=:branching[, t0])

Run a single multilevel simulation using the Gillespie algorithm, with the 
`population` giving the initial state. Simulation runs until the module population 
size reaches `maxmodules` or the age of the population reaches `tmax`. 

"""
function simulate_condfixtime!(population, tmax, maxmodules, birthrate, deathrate, moranrate, asymmetricrate, branchrate, 
    modulesize, branchinitsize, modulebranching, μ, mutationdist, moranincludeself, nextID, nextmoduleID, rng; moduleupdate=:branching, t0=nothing)

    t = isnothing(t0) ? age(population) : t0
    mosaic_mutations = Vector{Tuple{Int64, Float64}}[[]]
    condfixtimes = Vector{Float64}[[]]
    while t < tmax && (moduleupdate==:moran || length(population) < maxmodules)
        population, t, nextID, nextmoduleID = 
            update_population_condfixtimes!(population, mosaic_mutations, condfixtimes, birthrate, deathrate, moranrate, asymmetricrate, branchrate,  modulesize, 
                branchinitsize, modulebranching, t, nextID, nextmoduleID, μ, mutationdist, tmax, maxmodules, moranincludeself, rng; moduleupdate)
        #returns empty list of modules if population dies out
        if length(population) == 0
            return population, nextID, nextmoduleID
        end
    end
    return population, condfixtimes
end

function simulate_condfixtime_to_nfixed!(population, nfixed, tmax, maxmodules, birthrate, deathrate, moranrate, asymmetricrate, branchrate, 
    modulesize, branchinitsize, modulebranching, μ, mutationdist, moranincludeself, nextID, nextmoduleID, rng; moduleupdate=:branching, t0=nothing)
    t = isnothing(t0) ? age(population) : t0
    mosaic_mutations = Vector{Tuple{Int64, Float64}}[[]]
    condfixtimes = Vector{Float64}[[]]
    while t < tmax && (moduleupdate==:moran || length(population) < maxmodules)
        population, t, nextID, nextmoduleID = 
            update_population_condfixtimes!(population, mosaic_mutations, condfixtimes, birthrate, deathrate, moranrate, asymmetricrate, branchrate,  modulesize, 
                branchinitsize, modulebranching, t, nextID, nextmoduleID, μ, mutationdist, tmax, maxmodules, moranincludeself, rng; moduleupdate)
        #returns empty list of modules if population dies out
        if length(population) == 0
            return population, nextID, nextmoduleID
        end
        if length(condfixtimes[1]) >= nfixed
            return population, condfixtimes
        end
    end
    return population, condfixtimes
end

"""
    update_population!(population, birthrate, deathrate, moranrate, branchrate, modulesize, 
        branchinitsize, nextID, t, tmax, rng)

Perform a single Gillespie step to advance the simulation of a multilevel population, and
increase the time t → t + Δt.

Calculates the transition rates for a Moran, birth, death or module branching event; then 
selects from these transitions with probability proportional to rate. Time is increased by
`Δt ~ Exp(sum(transitionrates))`. If new time exceeds `tmax` do not perform any transition.
"""

function update_population_condfixtimes!(population, mosaic_mutations, condfixtimes, birthrate, deathrate, moranrate, asymmetricrate, branchrate, modulesize, 
    branchinitsize, modulebranching, t, nextID, nextmoduleID, μ, mutationdist, tmax, maxmodules, moranincludeself, rng; moduleupdate=:branching)
    
    transitionrates = get_transitionrates(population, birthrate, deathrate, 
        moranrate, asymmetricrate, branchrate, modulesize)
    t += exptime(rng, sum(transitionrates))
    #only update the population if t < tmax
    if t < tmax
        #choose transition type: 1=moran, 2=asymmetric, 3=birth, 4=death, 5=branch
        transitionid = sample(
            rng, 
            1:5, 
            ProbabilityWeights(transitionrates ./ sum(transitionrates))
        )
        population, nextID, nextmoduleID = transition_condfixtimes!(
            population, 
            mosaic_mutations, 
            condfixtimes,
            transitionid, 
            modulesize, 
            branchinitsize, 
            modulebranching, 
            t, 
            nextID,
            nextmoduleID,
            μ, 
            mutationdist,
            maxmodules,
            moranincludeself,
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
function transition_condfixtimes!(population, mosaic_mutations, condfixtimes, transitionid, modulesize, branchinitsize, 
    modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
    moranincludeself, rng; moduleupdate=:branching)
    
    if transitionid == 1
        _, nextID = moranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself)
    elseif transitionid == 2
        _, nextID = asymmetricupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 3
        _, nextID = birthupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 4
        deathupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, μ, mutationdist, rng)
    elseif transitionid == 5
        if moduleupdate == :branching || length(population) < maxmodules
            _, nextmoduleID = modulebranchingupdate_condfixtimes!(
                population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, branchinitsize, t, rng; 
                modulebranching, nextID, μ, mutationdist
            )
        elseif moduleupdate == :moran
            _, nextmoduleID = modulemoranupdate_condfixtimes!(
                population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, branchinitsize, t, rng; 
                modulebranching, nextID, μ, mutationdist
            )
        end
    end
    return population, nextID, nextmoduleID
end

function checkfixationextinction_addcondfixtimes!(mosaic_mutations_in_cellmodule, condfixtimes_in_cellmodule, cellmodule, t)
    ids_to_remove = Int64[]
    for (i, (mut, t0)) in enumerate(mosaic_mutations_in_cellmodule)
        fixed, extinct = isfixed_isextinct(mut, cellmodule)
        if fixed
            append!(condfixtimes_in_cellmodule, t - t0)
            push!(ids_to_remove, i)
        elseif extinct
            push!(ids_to_remove, i)
        end
    end
    deleteat!(mosaic_mutations_in_cellmodule, ids_to_remove)
end



function isfixed_isextinct(mutationid, cellmodule)
    mut_in_cell = [
        (mutationid in cell.mutations)
            for cell in cellmodule.cells
    ]
    return (all(mut_in_cell), !any(mut_in_cell))
end            

function new_mosaic_mutations!(mosaic_mutations_in_cellmodule, firstID, lastID, t)
    for mutID in firstID:lastID
        push!(mosaic_mutations_in_cellmodule, (mutID, t))
    end
end

"""
    moranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself=true)

Selects a homeostatic module uniformly at random to undergo a single Moran update. From 
that module one cell divides and one cell dies.
"""
function moranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself=true)
    cellmodule, parentcell, deadcell, cellmodule_id = 
        choose_homeostaticmodule_cells(population, rng; moranincludeself, maxmodulesize=modulesize)
    nextID_before = nextID
    _, nextID = celldivision!(cellmodule, parentcell, t, nextID, μ, mutationdist, rng)
    new_mosaic_mutations!(mosaic_mutations[cellmodule_id], nextID_before, nextID-1, t)
    celldeath!(cellmodule, deadcell, t, μ, mutationdist, rng)
    updatetime!(cellmodule, t)
    checkfixationextinction_addcondfixtimes!(mosaic_mutations[cellmodule_id], condfixtimes[cellmodule_id], cellmodule, t)
    return population, nextID
end

"""
    asymmetricupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)

Selects a homeostatic module uniformly at random to undergo a single asymmetric update. From 
that module one cell divides, producing a single offspring.
"""
function asymmetricupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)
    cellmodule, parentcell, _, cellmodule_id =
        choose_homeostaticmodule_cells(population, modulesize, rng; twocells=false, maxmodulesize=modulesize)
    nextID_before = nextID
    _, nextID = celldivision!(cellmodule, parentcell, t, nextID, μ, mutationdist, rng; nchildcells=1)
    new_mosaic_mutations!(mosaic_mutations[cellmodule_id], nextID_before, nextID-1, t)
    updatetime!(cellmodule, t)
    return population, nextID
end


"""
    birthupdate_condfixtimes!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)
    cellmodule, parentcell, cellmodule_id = 
        choose_growingmodule_cell(population, modulesize, rng)
    nextID_before = nextID
    _, nextID = celldivision!(cellmodule, parentcell, t, nextID, μ, mutationdist, rng)
    new_mosaic_mutations!(mosaic_mutations[cellmodule_id], nextID_before, nextID-1, t)
    updatetime!(cellmodule, t)
    return population, nextID
end


"""
    deathupdate!(population, modulesize, t, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to die. If 
cell death results in an empty module, remove that module from the population.
"""

function deathupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, μ, mutationdist, rng)
    cellmodule, deadcell, cellmodule_id = 
        choose_growingmodule_cell(population, modulesize, rng)
    celldeath!(cellmodule, deadcell, t, μ, mutationdist, rng)
    updatetime!(cellmodule, t)
    if length(cellmodule) == 0
        moduledeath!(population, cellmodule, t, μ, mutationdist, rng)
        deleteat!(mosaic_mutations, cellmodule)
    else
        checkfixation_addcondfixtimes!(mosaic_mutations[cellmodule_id], condfixtimes[cellmodule_id], cellmodule, t)
    end
    return population
end

"""
    modulebranchingupdate!(population, modulesize, branchinitsize, t, rng)
Select a homeostatic module, uniformly at random, to undergo branching. Cells are sampled 
(number of cells giving by `branchinitsize`) from the parent module to form a new module, 
which is added to `population`
"""
function modulebranchingupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, branchinitsize, t, rng; 
    modulebranching=:split, nextID=nothing, μ=nothing, mutationdist=nothing)
    
    parentmodule, parentmodule_id = choose_homeostaticmodule(population, modulesize, rng)
    parentmodule, newmodule, nextID = 
        newmoduleformation!(parentmodule, nextmoduleID, branchinitsize, t, rng; 
            modulebranching, nextID, μ, mutationdist)
    push!(population, newmodule)
    push!(mosaic_mutations, copy(mosaic_mutations[parentmodule_id]))
    push!(condfixtimes, copy(condfixtimes[parentmodule_id]))
    checkfixationextinction_addcondfixtimes!(mosaic_mutations[parentmodule_id], condfixtimes[parentmodule_id], parentmodule, t)
    checkfixationextinction_addcondfixtimes!(mosaic_mutations[end], condfixtimes[end], population[end], t)
    return population, nextmoduleID + 1
end

function modulemoranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, branchinitsize, t, rng; 
    modulebranching=:split, nextID=nothing, μ=nothing, mutationdist=nothing)

    parentmodule, parentmodule_id = choose_homeostaticmodule(population, modulesize, rng)
    parentmodule, newmodule, nextID = 
        newmoduleformation!(parentmodule, nextmoduleID, branchinitsize, t, rng; 
            modulebranching, nextID, μ, mutationdist)
    push!(population, newmodule)
    push!(mosaic_mutations, copy(mosaic_mutations[parentmodule_id]))
    push!(condfixtimes, copy(condfixtimes[parentmodule_id]))
    checkfixationextinction_addcondfixtimes!(mosaic_mutations[parentmodule_id], condfixtimes[parentmodule_id], parentmodule, t)
    checkfixationextinction_addcondfixtimes!(mosaic_mutations[end], condfixtimes[end], population[end], t)
    deadmodule_id = rand(rng, 1:length(population))
    deadmodule = population[deadmodule_id]
    moduledeath!(population, deadmodule, t, μ, mutationdist, rng)
    deleteat!(mosaic_mutations, deadmodule_id)
    deleteat!(condfixtimes, deadmodule_id)
    return population, nextmoduleID + 1
end