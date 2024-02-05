"""
    simulate_condfixtime!(population, tmax, maxmodules, birthrate, deathrate, moranrate, branchrate, 
        modulesize, branchinitsize, rng; moduleupdate=:branching[, t0])

Run a single multilevel simulation using the Gillespie algorithm, with the 
`population` giving the initial state. Simulation runs until the module population 
size reaches `maxmodules` or the age of the population reaches `tmax`. 

"""
function simulate_condfixtime!(population, input, ::NeutralSelection, counters, rng; timefunc)

    tmax = input.tmax
    nextID, nextmoduleID = counters
    t = age(population)
    transitionrates = get_neutral_transitionrates(population, input.branchrate, input.modulesize, input.quiescence)
    moduleupdate = getmoduleupdate(input)
    mosaic_mutations = (homeostatic=Vector{Tuple{Int64, Float64}}[[]],growing=Vector{Tuple{Int64, Float64}}[[]])
    condfixtimes = (homeostatic=Vector{Float64}[], growing=Vector{Float64}[[]])
    while t < tmax && (moduleupdate==:moran || length(population) < maxmodules)
        population, transitionrates, t, nextID, nextmoduleID = 
            update_population_neutral_condfixtimes!(
                population, 
                transitionrates,
                mosaic_mutations, 
                condfixtimes, 
                input.branchrate,  
                input.modulesize, 
                input.branchinitsize, 
                input.modulebranching, 
                input.quiescence,
                t, 
                nextID, 
                nextmoduleID, 
                input.μ, 
                input.mutationdist, 
                tmax, 
                input.maxmodules, 
                input.moranincludeself, 
                rng; 
                moduleupdate,
                timefunc)
        #returns empty list of modules if population dies out
        if length(population) == 0
            return population, nextID, nextmoduleID
        end
    end
    return population, condfixtimes
end

"""
    update_population_neutral_condfixtimes!(population, mosaic_mutations, condfixtimes, 
        branchrate, modulesize, branchinitsize, modulebranching, quiescence, t, nextID, 
        nextmoduleID, μ, mutationdist, tmax, maxmodules, moranincludeself, rng; 
        moduleupdate=:branching, timefunc=exptime)
"""
function update_population_neutral_condfixtimes!(population, transitionrates, mosaic_mutations, condfixtimes, 
    branchrate, modulesize, branchinitsize, modulebranching, quiescence, t, nextID, 
    nextmoduleID, μ, mutationdist, tmax, maxmodules, moranincludeself, rng; 
    moduleupdate=:branching, timefunc=exptime)
    
    t += exptime(rng, sum(transitionrates))
    #only update the population if t < tmax
    if t < tmax
        #choose transition type: 1=moran, 2=asymmetric, 3=birth, 4=death, 5=branch
        transitionid = sample(
            rng, 
            1:length(transitionrates), 
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
        if transitionid > 2
            update_neutral_transitionrates!(transitionrates, population, branchrate, modulesize)
        end
    end
    return population, transitionrates, t, nextID, nextmoduleID
end

"""
    transition_condfixtimes!(population, mosaic_mutations, condfixtimes, transitionid, modulesize, branchinitsize, 
        modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
        moranincludeself, rng; moduleupdate=:branching)

"""
function transition_condfixtimes!(population, mosaic_mutations, condfixtimes, transitionid, modulesize, branchinitsize, 
    modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
    moranincludeself, rng; moduleupdate=:branching)
    
    if transitionid == 1
        _, nextID = moranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself)
    elseif transitionid == 2
        _, nextID = asymmetricupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 3
        _, nextID = birthupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 4
        deathupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, μ, mutationdist, rng)
    elseif transitionid == 5
        if moduleupdate == :branching || length(population) < maxmodules
            _, nextmoduleID, nextID = modulebranchingupdate_condfixtimes!(
                population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, branchinitsize, t, rng; 
                modulebranching, nextID, μ, mutationdist
            )
        elseif moduleupdate == :moran
            _, nextmoduleID, nextID = modulemoranupdate_condfixtimes!(
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
"""
function moranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself=true)
    homeostaticmoduleid, parentcellid, deadcellid = 
        choose_homeostaticmodule_cells(population, rng; moranincludeself, maxmodulesize=modulesize)
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]    
    nextID_before = nextID
    _, _, nextID = celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, μ, mutationdist, rng)
    new_mosaic_mutations!(mosaic_mutations.homeostatic[homeostaticmoduleid], nextID_before, nextID-1, t)
    celldeath!(homeostaticmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
    updatetime!(homeostaticmodule, t)
    checkfixationextinction_addcondfixtimes!(
        mosaic_mutations.homeostatic[homeostaticmoduleid], 
        condfixtimes.homeostatic[homeostaticmoduleid], 
        homeostaticmodule, 
        t
    )
    return population, nextID
end

"""
    asymmetricupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)

Selects a homeostatic module uniformly at random to undergo a single asymmetric update. From 
that module one cell divides, producing a single offspring.
"""
function asymmetricupdate_condfixtimes!(population, mosaic_mutations, modulesize, t, nextID, μ, mutationdist, rng)
    homeostaticmoduleid, parentcellid = 
        choose_homeostaticmodule_cells(population, rng; twocells=false, maxmodulesize=modulesize)
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]    
    nextID_before = nextID
    _, _, nextID = celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, μ, mutationdist, rng; nchildcells=1)
    new_mosaic_mutations!(mosaic_mutations.homeostatic[homeostaticmoduleid], nextID_before, nextID-1, t)
    updatetime!(homeostaticmodule, t)
    return population, nextID
end


"""
    birthupdate_condfixtimes!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, maxmodulesize, t, nextID, μ, mutationdist, rng)
    growingmoduleid, parentcellid = choose_growingmodule_cell(population, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    nextID_before = nextID
    _, _, nextID = celldivision!(growingmodule, population.subclones, parentcellid, t, nextID, μ, mutationdist, rng)
    new_mosaic_mutations!(mosaic_mutations.growing[growingmoduleid], nextID_before, nextID-1, t)
    updatetime!(growingmodule, t)
    if length(growingmodule) >= maxmodulesize
        move_module_to_homeostasis!(population, growingmoduleid)
        push!(mosaic_mutations.homeostatic, popat!(mosaic_mutations.growing, growingmoduleid))
        push!(condfixtimes.homeostatic, popat!(condfixtimes.growing, growingmoduleid))
    end
    return population, nextID
end


"""
    deathupdate!(population, modulesize, t, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to die. If 
cell death results in an empty module, remove that module from the population.
"""

function deathupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, modulesize, t, μ, mutationdist, rng)
    growingmoduleid, deadcellid = choose_growingmodule_cell(population, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    celldeath!(growingmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
    updatetime!(growingmodule, t)
    if length(growingmodule) == 0
        moduledeath!(population, growingmodule, t, μ, mutationdist, rng)
        deleteat!(mosaic_mutations.growing, growingmoduleid)
        deleteat!(condfixtimes.growing, growingmoduleid)
    else
        checkfixation_addcondfixtimes!(mosaic_mutations.growing[growingmoduleid], condfixtimes.growing[growingmoduleid], growingmodule, t)
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
    
    parentmoduleid = choose_homeostaticmodule(population, rng)
    parentmodule = population.homeostatic_modules[parentmoduleid]
    parentmodule, newmodule, nextID = 
        newmoduleformation!(parentmodule, population.subclones, nextmoduleID, branchinitsize, t, rng; 
            modulebranching, nextID, μ, mutationdist)
    push!(population.growing_modules, newmodule)
    push!(mosaic_mutations.growing, copy(mosaic_mutations.homeostatic[parentmoduleid]))
    push!(condfixtimes.growing, copy(condfixtimes.homeostatic[parentmoduleid]))
    checkfixationextinction_addcondfixtimes!(mosaic_mutations.homeostatic[parentmoduleid], condfixtimes.homeostatic[parentmoduleid], parentmodule, t)
    checkfixationextinction_addcondfixtimes!(mosaic_mutations.growing[end], condfixtimes.growing[end], population[end], t)
    if modulebranching == :split
        move_module_a_to_b!(population.homeostatic_modules, population.growing_modules, parentmoduleid)
        push!(mosaic_mutations.growing, popat!(mosaic_mutations.homeostatic, parentmoduleid))
        push!(condfixtimes.growing, popat!(condfixtimes.homeostatic, parentmoduleid))
    end
    return population, nextmoduleID + 1, nextID
end

function modulemoranupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, branchinitsize, t, rng; 
    modulebranching=:split, nextID=nothing, μ=nothing, mutationdist=nothing)

    _, _, nextID = modulebranchingupdate_condfixtimes!(population, mosaic_mutations, condfixtimes, nextmoduleID, modulesize, 
        branchinitsize, t, rng; modulebranching, nextID, μ, mutationdist)

    deadmoduleid = choose_any_module(population, rng)
    deadmodule = population[deadmoduleid]
    Nhom = length(population.homeostatic_modules)
    moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
    if deadmoduleid <= Nhom
        deleteat!(mosaic_mutations.homeostatic, deadmoduleid)
        deleteat!(condfixtimes.homeostatic, deadmoduleid)
    else
        deleteat!(mosaic_mutations.growing, deadmoduleid - Nhom)
        deleteat!(condfixtimes.growing, deadmoduleid - Nhom)
    end
    return population, nextmoduleID + 1, nextID
end