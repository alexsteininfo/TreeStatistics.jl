###
### Implementation for multilevel neutral simulations
###

"""
    simulate!(population, input, selection::NeutralSelection, counters, rng;
    timefunc=exptime, t0=nothing, tmax=Inf)

Run a single multilevel simulation using the Gillespie algorithm, with the
`population` giving the initial state. Simulation runs until the module population
size reaches `maxmodules` or the age of the population reaches `tmax`.

"""
function simulate!(
    population,
    input::MultilevelInput,
    ::NeutralSelection,
    counters,
    rng;
    timefunc=exptime,
    t0=nothing,
    tmax=Inf
)
    tmax = minimum((input.tmax, tmax))
    nextID, nextmoduleID = counters
    t = isnothing(t0) ? age(population) : t0
    seasonalstate = initializeseason(input.quiescence)
    transitionrates = get_neutral_transitionrates(
        population,
        input.branchrate,
        input.modulesize,
        input.quiescence,
        seasonalstate
    )
    moduleupdate = getmoduleupdate(input)
    while t < tmax && (moduleupdate==:moran || length(population) < input.maxmodules)
        population, transitionrates, t, nextID, nextmoduleID =
            update_population_neutral!(
                population,
                transitionrates,
                input.branchrate,
                input.modulesize,
                input.branchinitsize,
                input.modulebranching,
                input.quiescence,
                seasonalstate,
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
                timefunc
            )
        #returns empty list of modules if population dies out
        if length(population) == 0
            return population, nextID, nextmoduleID
        end
    end
    #add final time-dependent mutations
    final_timedep_mutations!(
        population,
        input.μ,
        input.mutationdist,
        rng;
        tend=tmax,
        mutID=nextID
    )
    return population, nextID, nextmoduleID
end

getmoduleupdate(::MultilevelBranchingInput) = :branching
getmoduleupdate(::MultilevelBranchingMoranInput) = :moran
getmoduleupdate(::MultilevelMoranInput) = :moran

"""
    update_population_neutral!(
        population,
        transitionrates,
        branchrate,
        modulesize,
        branchinitsize,
        modulebranching,
        quiescence,
        seasonalstate,
        t,
        nextID,
        nextmoduleID,
        μ,
        mutationdist,
        tmax,
        maxmodules,
        moranincludeself,
        rng;
        moduleupdate=:branching,
        timefunc=exptime
    )

Perform a single Gillespie step to advance the simulation of a multilevel population, and
increase the time t → t + Δt.

Calculates the transition rates for a Moran, birth, death or module branching event; then
selects from these transitions with probability proportional to rate. Time is increased by
`Δt ~ Exp(sum(transitionrates))`. If new time exceeds `tmax` do not perform any transition.
"""
function update_population_neutral!(
    population,
    transitionrates,
    branchrate,
    modulesize,
    branchinitsize,
    modulebranching,
    quiescence,
    seasonalstate,
    t,
    nextID,
    nextmoduleID,
    μ,
    mutationdist,
    tmax,
    maxmodules,
    moranincludeself,
    rng;
    moduleupdate=:branching,
    timefunc=exptime
)
    t += timefunc(rng, sum(transitionrates))
    #only update the population if t < tmax
    if t < tmax
        #choose transition type: 1=moran, 2=asymmetric, 3=birth, 4=death, 5=branch
        transitionid = sample(
            rng,
            1:length(transitionrates),
            ProbabilityWeights(transitionrates ./ sum(transitionrates))
        )
        population, nextID, nextmoduleID = transition!(
            population,
            transitionid,
            modulesize,
            branchinitsize,
            modulebranching,
            quiescence,
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
        seasonswitch, seasonalstate = switchseasons(t, quiescence, seasonalstate)
        if transitionupdaterequired(transitionid, quiescence, seasonswitch)
            update_neutral_transitionrates!(
                transitionrates,
                population,
                branchrate,
                modulesize,
                quiescence,
                seasonalstate
            )
        end
    end
    return population, transitionrates, t, nextID, nextmoduleID
end

initializeseason(::Any) = nothing
function initializeseason(quiescence::SeasonalQuiescence)
    return (winter = false, nextswitchtime = quiescence.summerduration)
end

switchseasons(t, ::NoQuiescence, ::Any) = (false, nothing)
switchseasons(t, ::StochasticQuiescence, ::Any) = (false, nothing)
function switchseasons(t, quiescence::SeasonalQuiescence, seasonalstate)
    if t > seasonalstate.nextswitchtime
        seasonalstate =
            if seasonalstate.winter
                nextswitchtime = seasonalstate.nextswitchtime + quiescence.summerduration
                (winter = false, nextswitchtime = nextswitchtime)
            else
                nextswitchtime = seasonalstate.nextswitchtime + quiescence.winterduration
                (winter = true, nextswitchtime = nextswitchtime)
            end
        return (true, seasonalstate)
    else
        return (false, seasonalstate)
    end
end

transitionupdaterequired(transitionid, ::NoQuiescence, ::Any) = transitionid > 2
function transitionupdaterequired(transitionid, ::SeasonalQuiescence, seasonswitch)
    transitionid > 2 || seasonswitch
end
transitionupdaterequired(transitionid, ::StochasticQuiescence, ::Any) = transitionid > 4

"""
    transition!(
        population,
        transitionid,
        modulesize,
        branchinitsize,
        modulebranching,
        ::AbstractDeterministicQuiescence,
        t,
        nextID,
        nextmoduleID,
        μ,
        mutationdist,
        maxmodules,
        moranincludeself,
        rng;
        moduleupdate=:branching
    )

Perform a single transition step on `population`, determined by `transitionid`.
"""
function transition!(
    population,
    transitionid,
    modulesize,
    branchinitsize,
    modulebranching,
    ::AbstractDeterministicQuiescence,
    t,
    nextID,
    nextmoduleID,
    μ,
    mutationdist,
    maxmodules,
    moranincludeself,
    rng;
    moduleupdate=:branching
)
    if transitionid == 1
        _, nextID = moranupdate!(
            population,
            modulesize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            moranincludeself
        )
    elseif transitionid == 2
        _, nextID = asymmetricupdate!(
            population,
            modulesize,
            t,
            nextID,
            μ,
            mutationdist,
            rng
        )
    elseif transitionid == 3
        _, nextID = birthupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 4
        deathupdate!(population, t, μ, mutationdist, rng)
    elseif transitionid == 5
        if moduleupdate == :branching || length(population) < maxmodules
            _, nextmoduleID, nextID = modulebranchingupdate!(
                population, nextmoduleID, branchinitsize, t, rng;
                modulebranching, nextID, μ, mutationdist
            )
        elseif moduleupdate == :moran
            _, nextmoduleID, nextID = modulemoranupdate!(
                population, nextmoduleID, branchinitsize, t, rng;
                modulebranching, nextID, μ, mutationdist
            )
        end
    end
    return population, nextID, nextmoduleID
end

function transition!(
    population::PopulationWithQuiescence,
    transitionid,
    modulesize,
    branchinitsize,
    modulebranching,
    ::StochasticQuiescence,
    t,
    nextID,
    nextmoduleID,
    μ,
    mutationdist,
    maxmodules,
    moranincludeself,
    rng;
    moduleupdate=:branching
)

    if transitionid == 1
        _, nextID = moranupdate!(
            population,
            modulesize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            moranincludeself
        )
    elseif transitionid == 2
        _, nextID = moranupdate_quiescent!(
            population,
            modulesize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            moranincludeself
        )
    elseif transitionid == 3
        _, nextID = asymmetricupdate!(
            population,
            modulesize,
            t,
            nextID,
            μ,
            mutationdist,
            rng
        )
    elseif transitionid == 4
        _, nextID = asymmetricupdate_quiescent!(
            population,
            modulesize,
            t,
            nextID,
            μ,
            mutationdist,
            rng
        )
    elseif transitionid == 5
        _, nextID = birthupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    elseif transitionid == 6
        deathupdate!(population, t, μ, mutationdist, rng)
    elseif transitionid == 7
        if moduleupdate == :branching || length(population) < maxmodules
            _, nextmoduleID, nextID = modulebranchingupdate!(
                population, nextmoduleID, branchinitsize, t, rng;
                modulebranching, nextID, μ, mutationdist
            )
        elseif moduleupdate == :moran
            _, nextmoduleID, nextID = modulemoranupdate!(
                population, nextmoduleID, branchinitsize, t, rng;
                modulebranching, nextID, μ, mutationdist
            )
        end
    elseif transitionid == 8
        if moduleupdate == :branching || length(population) < maxmodules
            _, nextmoduleID, nextID = modulebranchingupdate_quiescent!(
                population, nextmoduleID, branchinitsize, t, rng;
                modulebranching, nextID, μ, mutationdist
            )
        elseif moduleupdate == :moran
            _, nextmoduleID, nextID = modulemoranupdate_quiescent!(
                population, nextmoduleID, branchinitsize, t, rng;
                modulebranching, nextID, μ, mutationdist
            )
        end
    elseif transitionid == 9
        quiescenceonupdate!(population, rng)
    elseif transitionid == 10
        quiescenceoffupdate!(population, rng)
    end
    return population, nextID, nextmoduleID
end

"""
    moranupdate!(
        population,
        modulesize,
        t,
        nextID,
        μ,
        mutationdist,
        rng;
        moranincludeself=true
)

Selects a homeostatic module uniformly at random to undergo a single Moran update. From
that module one cell divides and one cell dies.
"""
function moranupdate!(
    population,
    modulesize,
    t,
    nextID,
    μ,
    mutationdist,
    rng;
    moranincludeself=true
)
    homeostaticmoduleid, parentcellid, deadcellid =
        choose_homeostaticmodule_cells(
            population,
            rng;
            moranincludeself,
            maxmodulesize=modulesize
        )
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
    nextID = _moranupdate!(
        homeostaticmodule,
        population.subclones,
        parentcellid,
        deadcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng
    )

    return population, nextID
end

function moranupdate_quiescent!(
    population,
    modulesize,
    t,
    nextID,
    μ,
    mutationdist,
    rng;
    moranincludeself=true
)
    quiescentmoduleid, parentcellid, deadcellid =
        choose_quiescentmodule_cells(
            population,
            rng;
            moranincludeself,
            maxmodulesize=modulesize
        )
    quiescentmodule = population.quiescent_modules[quiescentmoduleid]
    nextID = _moranupdate!(
        quiescentmodule,
        population.subclones,
        parentcellid,
        deadcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng
    )

    return population, nextID
end

function _moranupdate!(
    chosenmodule,
    subclones,
    parentcellid,
    deadcellid,
    t,
    nextID,
    μ,
    mutationdist,
    rng
)
    _, _, nextID = celldivision!(
        chosenmodule,
        subclones,
        parentcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng
    )
    celldeath!(chosenmodule, subclones, deadcellid, t, μ, mutationdist, rng)
    updatetime!(chosenmodule, t)
    return nextID
end

"""
    asymmetricupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a homeostatic module uniformly at random to undergo a single asymmetric update. From
that module one cell divides, producing a single offspring.
"""
function asymmetricupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    homeostaticmoduleid, parentcellid =
        choose_homeostaticmodule_cells(
            population,
            rng;
            twocells=false,
            maxmodulesize=modulesize
        )
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
    nextID = _asymmetricupdate!(
        homeostaticmodule,
        population.subclones,
        parentcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng
    )
    return population, nextID
end

function asymmetricupdate_quiescent!(
    population,
    modulesize,
    t,
    nextID,
    μ,
    mutationdist,
    rng
)
    quiescentmoduleid, parentcellid = choose_quiescentmodule_cells(
        population,
        rng;
        twocells=false,
        maxmodulesize=modulesize
    )
    quiescentmodule = population.quiescent_modules[quiescentmoduleid]
    nextID = _asymmetricupdate!(
        quiescentmodule,
        population.subclones,
        parentcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng
    )
    return population, nextID
end

function _asymmetricupdate!(
    chosenmodule,
    subclones,
    parentcellid,
    t,
    nextID,
    μ,
    mutationdist,
    rng
)
    _, _, nextID = celldivision!(
        chosenmodule,
        subclones,
        parentcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng;
        nchildcells=1
    )
    updatetime!(chosenmodule, t)
    return nextID
end

"""
    birthupdate!(population, maxmodulesize, t, nextID, μ, mutationdist, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate!(population, maxmodulesize, t, nextID, μ, mutationdist, rng)
    growingmoduleid, parentcellid = choose_growingmodule_cell(population, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    _, _, nextID = celldivision!(
        growingmodule,
        population.subclones,
        parentcellid,
        t,
        nextID,
        μ,
        mutationdist,
        rng
    )
    updatetime!(growingmodule, t)
    if length(growingmodule) == maxmodulesize
        move_module_to_homeostasis!(population, growingmoduleid)
    end
    return population, nextID
end

"""
    deathupdate!(population, modulesize, t, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to die. If
cell death results in an empty module, remove that module from the population.
"""

function deathupdate!(population, t, μ, mutationdist, rng)
    growingmoduleid, deadcellid = choose_growingmodule_cell(population, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    celldeath!(growingmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
    updatetime!(growingmodule, t)
    if length(growingmodule) == 0
        moduledeath!(population, growingmoduleid; moduletype=:growing)
    end
    return population
end

"""
    modulebranchingupdate!(population, modulesize, branchinitsize, t, rng)
Select a homeostatic module, uniformly at random, to undergo branching. Cells are sampled
(number of cells giving by `branchinitsize`) from the parent module to form a new module,
which is added to `population`
"""
function modulebranchingupdate!(population, nextmoduleID, branchinitsize, t, rng;
    modulebranching=:split, nextID=nothing, μ=nothing, mutationdist=nothing)

    parentmoduleid = choose_homeostaticmodule(population, rng)
    nextID =
        _modulebranchingupdate!(
            parentmoduleid,
            population.homeostatic_modules,
            population,
            nextmoduleID,
            branchinitsize,
            t,
            rng,
            modulebranching,
            nextID,
            μ,
            mutationdist
        )
    return population, nextmoduleID + 1, nextID
end

function modulebranchingupdate_quiescent!(
    population,
    nextmoduleID,
    branchinitsize,
    t,
    rng;
    modulebranching=:split,
    nextID=nothing,
    μ=nothing,
    mutationdist=nothing
)
    parentmoduleid = choose_quiescentmodule(population, rng)
    nextID =
        _modulebranchingupdate!(
            parentmoduleid,
            population.quiescent_modules,
            population,
            nextmoduleID,
            branchinitsize,
            t,
            rng,
            modulebranching,
            nextID,
            μ,
            mutationdist
        )
    return population, nextmoduleID + 1, nextID
end

function _modulebranchingupdate!(
    parentmoduleid,
    modulelist,
    population,
    nextmoduleID,
    branchinitsize,
    t,
    rng,
    modulebranching,
    nextID,
    μ,
    mutationdist
)
    parentmodule = modulelist[parentmoduleid]
    parentmodule, newmodule, nextID = newmoduleformation!(
        parentmodule,
        population.subclones,
        nextmoduleID,
        branchinitsize,
        t,
        rng;
        modulebranching,
        nextID,
        μ,
        mutationdist
    )
    push!(population.growing_modules, newmodule)
    if modulebranching == :split
        move_module_a_to_b!(modulelist, population.growing_modules, parentmoduleid)
    end
    return nextID
end

function modulemoranupdate!(
    population,
    nextmoduleID,
    branchinitsize,
    t,
    rng;
    modulebranching=:split,
    nextID=nothing,
    μ=nothing,
    mutationdist=nothing
)
    _, _, nextID = modulebranchingupdate!(
        population,
        nextmoduleID,
        branchinitsize,
        t,
        rng;
        modulebranching,
        nextID,
        μ,
        mutationdist
    )
    deadmoduleid = choose_any_module(population, rng)
    moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
    return population, nextmoduleID + 1, nextID
end

function modulemoranupdate_quiescent!(
    population,
    nextmoduleID,
    branchinitsize,
    t,
    rng;
    modulebranching=:split,
    nextID=nothing,
    μ=nothing,
    mutationdist=nothing
)

    _, _, nextID = modulebranchingupdate_quiescent!(
        population,
        nextmoduleID,
        branchinitsize,
        t,
        rng;
        modulebranching,
        nextID,
        μ,
        mutationdist
    )
    deadmoduleid = choose_any_module(population, rng)
    moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
    return population, nextmoduleID + 1, nextID
end

function quiescenceonupdate!(population, rng)
    homeostaticmoduleid = choose_homeostaticmodule(population, rng)
    move_module_a_to_b!(
        population.homeostatic_modules,
        population.quiescent_modules,
        homeostaticmoduleid
    )
    return population
end

function quiescenceoffupdate!(population, rng)
    quiescentmoduleid = choose_quiescentmodule(population, rng)
    move_module_a_to_b!(
        population.quiescent_modules,
        population.homeostatic_modules,
        quiescentmoduleid
    )
    return population
end


function newmoduleformation!(
    parentmodule,
    subclones,
    nextmoduleID,
    branchinitsize::Int,
    t,
    rng;
    modulebranching=:split,
    nextID=nothing,
    μ=nothing,
    mutationdist=nothing
)
    if modulebranching == :withreplacement
        return sample_new_module_with_replacement!(
            parentmodule,
            subclones,
            nextmoduleID,
            branchinitsize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            timedepmutationsonly=false
        )
    elseif modulebranching == :withoutreplacement
        return sample_new_module_without_replacement!(
            parentmodule,
            subclones,
            nextmoduleID,
            branchinitsize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            timedepmutationsonly=false
        )
    elseif modulebranching == :withreplacement_nomutations
        return sample_new_module_with_replacement!(
            parentmodule,
            subclones,
            nextmoduleID,
            branchinitsize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            timedepmutationsonly=true
        )
    elseif modulebranching == :withoutreplacement_nomutations
        return sample_new_module_without_replacement!(
            parentmodule,
            subclones,
            nextmoduleID,
            branchinitsize,
            t,
            nextID,
            μ,
            mutationdist,
            rng;
            timedepmutationsonly=true
        )
    elseif modulebranching == :split
        cellmodule, newcellmodule =
            sample_new_module_split!(
                parentmodule,
                nextmoduleID,
                branchinitsize,
                t,
                rng
            )
        return cellmodule, newcellmodule, nextID
    else error("$modulebranching is not a valid module branching option")
    end
end

"""
    sample_new_module_with_replacement!(
        cellmodule,
        subclones,
        nextmoduleID,
        branchinitsize,
        branchtime,
        nextID,
        μ,
        mutationdist,
        rng::AbstractRNG;
        timedepmutationsonly=false
)

Sample cells with replacement, uniformly at random, from `cellmodule`. Each of these cells
divides, with one offspring remaining in `cellmodule`, the other becoming part of the new
module.

The number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.
"""
function sample_new_module_with_replacement!(
    cellmodule,
    subclones,
    nextmoduleID,
    branchinitsize,
    branchtime,
    nextID,
    μ,
    mutationdist,
    rng::AbstractRNG;
    timedepmutationsonly=false
)
    ncells = length(cellmodule.cells)
    newcells = eltype(cellmodule.cells)[]
    for i in 1:branchinitsize
        randcellid = rand(rng, 1:ncells)
        cellmodule, subclones, nextID = celldivision!(
            cellmodule,
            subclones,
            randcellid,
            branchtime,
            nextID,
            μ,
            mutationdist,
            rng;
            timedepmutationsonly
        )
        newcell = pop!(cellmodule.cells)
        push!(newcells, newcell)
    end

    newcellmodule = new_module_from_cells(
        newcells,
        branchtime,
        [branchtime],
        nextmoduleID,
        cellmodule.id,
        cellmodule.structure
        )
    updatetime!(cellmodule, branchtime)
    push!(cellmodule.branchtimes, branchtime)
    return cellmodule, newcellmodule, nextID
end

"""
    sample_new_module_without_replacement!(
        cellmodule,
        subclones,
        nextmoduleID,
        branchinitsize,
        branchtime,
        nextID,
        μ,
        mutationdist,
        rng::AbstractRNG;
        timedepmutationsonly=false
    )

Sample cells without replacement, uniformly at random, from `cellmodule`. Each of these
cells divides, with one offspring remaining in `cellmodule`, the other becoming part of the
new module.

The number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module_without_replacement!(
    cellmodule,
    subclones,
    nextmoduleID,
    branchinitsize,
    branchtime,
    nextID,
    μ,
    mutationdist,
    rng::AbstractRNG;
    timedepmutationsonly=false
)

    sampleids = sample(rng, 1:length(cellmodule.cells), branchinitsize, replace=false)
    newcells = eltype(cellmodule.cells)[]
    for cellid in sampleids
        cellmodule, subclones, nextID =
            celldivision!(cellmodule, subclones, cellid, branchtime, nextID, μ, mutationdist, rng;
                timedepmutationsonly)
        newcell = pop!(cellmodule.cells)
        push!(newcells, newcell)
    end

    newcellmodule =
        new_module_from_cells(newcells, branchtime, [branchtime], nextmoduleID,
            cellmodule.id, cellmodule.structure)
    updatetime!(cellmodule, branchtime)
    push!(cellmodule.branchtimes, branchtime)

    return cellmodule, newcellmodule, nextID
end

"""
    sample_new_module_split!(
        cellmodule,
        nextmoduleID,
        branchinitsize,
        branchtime,
        rng::AbstractRNG
    )

Sample cells without replacement, uniformly at random, from `cellmodule` to form the new
module. Remaining cells form the parent module.

The number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module_split!(
    cellmodule,
    nextmoduleID,
    branchinitsize,
    branchtime,
    rng::AbstractRNG
)

    sampleids = sample(rng, 1:length(cellmodule.cells), branchinitsize, replace=false)
    newcellmodule =
        new_module_from_cells(cellmodule.cells[sampleids], branchtime, [branchtime], nextmoduleID,
            cellmodule.id, cellmodule.structure)
    updatetime!(cellmodule, branchtime)
    push!(cellmodule.branchtimes, branchtime)
    cellremoval!(cellmodule, sampleids)

    return cellmodule, newcellmodule
end


"""
    choose_homeostaticmodule(population, rng::AbstractRNG)
Select a homeostatic module, i.e. module of size `maxmodulesize`, uniformly at random.
"""
function choose_homeostaticmodule(population, rng::AbstractRNG)
    homeostaticmoduleid = rand(rng, 1:length(population.homeostatic_modules))
    return homeostaticmoduleid
end

"""
    choose_quiescentmodule(population, rng::AbstractRNG)
Select a homeostatic module, i.e. module of size `maxmodulesize`, uniformly at random.
"""
function choose_quiescentmodule(population, rng::AbstractRNG)
    quiescentmoduleid = rand(rng, 1:length(population.quiescent_modules))
    return quiescentmoduleid
end

"""
    choose_any_module(population, rng::AbstractRNG)
Select a module uniformly at random from all modules. Returns `moduleid` and the vector it
is in (either `growing_modules` or `homeostatic_modules`.)
"""
function choose_any_module(population, rng::AbstractRNG)
    moduleid = rand(rng, 1:length(population))
    return moduleid
end

"""
    choose_homeostaticmodule_cells(
        population,
        rng::AbstractRNG;
        twocells=true,
        moranincludeself=true,
        maxmodulesize=length(first(population.homeostatic_modules))
    )
Select a homeostatic module and the ids of two cells from the module, uniformly at random.
"""
function choose_homeostaticmodule_cells(
    population,
    rng::AbstractRNG;
    twocells=true,
    moranincludeself=true,
    maxmodulesize=length(first(population.homeostatic_modules))
)
    chosenmoduleid = choose_homeostaticmodule(population, rng)
    dividecellid = rand(rng, 1:maxmodulesize)
    deadcellid =
        if !twocells
            nothing
        else
            choose_moran_deadcell(maxmodulesize, dividecellid, moranincludeself, rng)
        end
    return chosenmoduleid, dividecellid, deadcellid
end

"""
    choose_quiescentmodule_cells(
        population,
        rng::AbstractRNG;
        twocells=true,
        moranincludeself=true,
        maxmodulesize=length(first(population.quiescent_modules))
    )
Select a homeostatic module and the ids of two cells from the module, uniformly at random.
"""
function choose_quiescentmodule_cells(
    population,
    rng::AbstractRNG;
    twocells=true,
    moranincludeself=true,
    maxmodulesize=length(first(population.quiescent_modules))
)
    chosenmoduleid = choose_quiescentmodule(population, rng)
    dividecellid = rand(rng, 1:maxmodulesize)
    deadcellid =
        if !twocells
            nothing
        else
            choose_moran_deadcell(maxmodulesize, dividecellid, moranincludeself, rng)
        end
    return chosenmoduleid, dividecellid, deadcellid
end

"""
    choose_growingmodule_cell(population, rng::AbstractRNG)

Select a cell uniformly at random from all cells in growing (non-homeostatic) modules, and
return the module and cell id.
"""
function choose_growingmodule_cell(population, rng::AbstractRNG)
    growingmodulesizes = map(length, population.growing_modules)
    chosenmoduleid = sample(
        rng,
        1:length(growingmodulesizes),
        ProbabilityWeights(growingmodulesizes ./ sum(growingmodulesizes))
    )
    chosencellid = rand(rng, 1:growingmodulesizes[chosenmoduleid])
    return chosenmoduleid, chosencellid
end

"""
    get_neutral_transitionrates(
        population,
        branchrate,
        modulesize,
        quiescence=NoQuiescence(),
        seasonalstate=nothing
    )

Compute the rates for moran update, asymmetric update, birth update, death update and
module branching.
"""
function get_neutral_transitionrates(
    population,
    branchrate,
    modulesize,
    quiescence=NoQuiescence(),
    seasonalstate=nothing
)
    rates = zeros(Float64, number_neutral_transitions(quiescence))
    return update_neutral_transitionrates!(
        rates,
        population,
        branchrate,
        modulesize,
        quiescence,
        seasonalstate
    )
end

number_neutral_transitions(::AbstractDeterministicQuiescence) = 5
number_neutral_transitions(::StochasticQuiescence) = 10

function update_neutral_transitionrates!(rates, population, branchrate, modulesize,
    quiescence=NoQuiescence(), seasonalstate=nothing)

    cellrates = getwildtyperates(population)
    current_cellrates = getcurrentcellrates(cellrates, quiescence, seasonalstate)
    current_branchrate = getcurrentbranchrate(branchrate, quiescence, seasonalstate)
    number_homeostatic_modules = length(population.homeostatic_modules)
    number_cells_in_growing_modules = sum(length.(population.growing_modules))
    rates[1] = number_homeostatic_modules * current_cellrates.moranrate * modulesize
    rates[2] = number_homeostatic_modules * current_cellrates.asymmetricrate * modulesize
    rates[3] = number_cells_in_growing_modules * current_cellrates.birthrate
    rates[4] = number_cells_in_growing_modules * current_cellrates.deathrate
    rates[5] = number_homeostatic_modules * current_branchrate
    return rates
end

function update_neutral_transitionrates!(
    rates,
    population::PopulationWithQuiescence,
    branchrate,
    modulesize,
    quiescence::StochasticQuiescence,
    ::Any
)
    cellrates = getwildtyperates(population)
    number_homeostatic_modules = length(population.homeostatic_modules)
    number_quiescent_modules = length(population.quiescent_modules)
    number_cells_in_growing_modules = sum(length.(population.growing_modules))
    rates[1] = number_homeostatic_modules * cellrates.moranrate * modulesize
    rates[2] = (number_quiescent_modules * cellrates.moranrate
        * quiescence.divisionfactor * modulesize)
    rates[3] = number_homeostatic_modules * cellrates.asymmetricrate * modulesize
    rates[4] = (number_quiescent_modules * cellrates.asymmetricrate
       * quiescence.divisionfactor * modulesize)
    rates[5] = number_cells_in_growing_modules * cellrates.birthrate
    rates[6] = number_cells_in_growing_modules * cellrates.deathrate
    rates[7] = number_homeostatic_modules * branchrate
    rates[8] = number_quiescent_modules * branchrate * quiescence.branchfactor
    rates[9] = number_homeostatic_modules * quiescence.onrate
    rates[10] = number_quiescent_modules * quiescence.offrate
    return rates
end

getcurrentbranchrate(branchrate, ::NoQuiescence, ::Any) = branchrate

function getcurrentbranchrate(branchrate, quiescence::SeasonalQuiescence, seasonalstate)
    if seasonalstate.winter
        return branchrate * quiescence.branchfactor
    else
        return branchrate
    end
end

getcurrentcellrates(cellrates, ::NoQuiescence, ::Any) = cellrates

function getcurrentcellrates(cellrates, quiescence::SeasonalQuiescence, seasonalstate)
    if seasonalstate.winter
        return (
            birthrate = cellrates.birthrate * quiescence.divisionfactor,
            deathrate = cellrates.deathrate * quiescence.divisionfactor,
            moranrate = cellrates.moranrate * quiescence.divisionfactor,
            asymmetricrate = cellrates.asymmetricrate * quiescence.divisionfactor
        )
    else
        return cellrates
    end
end

# function addquiescencerates!(rates, quiescence::StochasticQuiescence)

# end

"""
    moduledeath!(population, cellmodule)
Kill all live cells in `cellmodule` and remove it from `population`.
    """
function moduledeath!(population, dyingmoduleid, t, μ, mutationdist, rng; moduletype=:all)
    population, dyingmodule = removemodule!(population, dyingmoduleid; moduletype)
    for cell in dyingmodule.cells
        population.subclones[getclonetype(cell)].size -= 1
    end
    length(dyingmodule) == 0 || killallcells!(dyingmodule.cells, t, μ, mutationdist, rng)
    return population
end

function moduledeath!(population, dyingmoduleid; moduletype=:all)
    population, dyingmodule = removemodule!(population, dyingmoduleid; moduletype)
    @assert length(dyingmodule) == 0 "trying to kill non-empty module"
    return population
end

function removemodule!(population::Population, dyingmoduleid; moduletype=:all)
    if moduletype == :homeostatic
        dyingmodule = population.homeostatic_modules[dyingmoduleid]
        deleteat!(population.homeostatic_modules, dyingmoduleid)
        return population, dyingmodule
    elseif moduletype == :growing
        dyingmodule = population.growing_modules[dyingmoduleid]
        deleteat!(population.growing_modules, dyingmoduleid)
        return population, dyingmodule
    elseif moduletype == :all
        Nhom = length(population.homeostatic_modules)
        if dyingmoduleid <= Nhom
            dyingmodule = population.homeostatic_modules[dyingmoduleid]
            deleteat!(population.homeostatic_modules, dyingmoduleid)
            return population, dyingmodule
        else
            dyingmodule = population.growing_modules[dyingmoduleid - Nhom]
            deleteat!(population.growing_modules, dyingmoduleid - Nhom)
            return population, dyingmodule
        end
    else error("$moduletype not an allowed `moduletype` option")
    end
end


function removemodule!(population::PopulationWithQuiescence, dyingmoduleid; moduletype=:all)
    if moduletype == :homeostatic
        dyingmodule = population.homeostatic_modules[dyingmoduleid]
        deleteat!(population.homeostatic_modules, dyingmoduleid)
        return population, dyingmodule
    elseif moduletype == :quiescent
        dyingmodule = population.quiescent_modules[dyingmoduleid]
        deleteat!(population.quiescent_modules, dyingmoduleid)
        return population, dyingmodule
    elseif moduletype == :growing
        dyingmodule = population.growing_modules[dyingmoduleid]
        deleteat!(population.growing_modules, dyingmoduleid)
        return population, dyingmodule
    elseif moduletype == :all
        Nhom = length(population.homeostatic_modules)
        Nqui = length(population.quiescent_modules)
        if dyingmoduleid <= Nhom
            dyingmodule = population.homeostatic_modules[dyingmoduleid]
            deleteat!(population.homeostatic_modules, dyingmoduleid)
            return population, dyingmodule
        elseif dyingmoduleid <= Nqui + Nhom
            dyingmodule = population.quiescent_modules[dyingmoduleid - Nhom]
            deleteat!(population.quiescent_modules, dyingmoduleid - Nhom)
            return population, dyingmodule
        else
            dyingmodule = population.growing_modules[dyingmoduleid - Nhom - Nqui]
            deleteat!(population.growing_modules, dyingmoduleid - Nhom - Nqui)
            return population, dyingmodule
        end
    else error("$moduletype not an allowed `moduletype` option")
    end
end

function move_module_to_homeostasis!(population, cellmoduleid::Integer)
    return move_module_a_to_b!(
        population.growing_modules, population.homeostatic_modules, cellmoduleid
    )
end

function move_module_to_growing!(population, cellmoduleid::Integer)
    return move_module_a_to_b!(
        population.homeostatic_modules, population.growing_modules, cellmoduleid
    )
end

function move_module_a_to_b!(modules_a, modules_b, cellmoduleid::Integer)
    cellmodule = popat!(modules_a, cellmoduleid)
    push!(modules_b, cellmodule)
    return modules_a, modules_b
end

killallcells!(population::CellVector, args...) = nothing
