"""
    simulate!(population, input, selection::NeutralSelection, counters, rng; 
    timefunc=exptime, t0=nothing, tmax=Inf)

Run a single multilevel simulation using the Gillespie algorithm, with the 
`population` giving the initial state. Simulation runs until the module population 
size reaches `maxmodules` or the age of the population reaches `tmax`. 

"""
function simulate!(population, input::MultilevelInput, ::NeutralSelection, counters, rng; 
    timefunc=exptime, t0=nothing, tmax=Inf)

    tmax = minimum((input.tmax, tmax))
    nextID, nextmoduleID = counters
    t = isnothing(t0) ? age(population) : t0
    transitionrates = get_neutral_transitionrates(population, input.branchrate, input.modulesize)
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
    return population, nextID, nextmoduleID
end

getmoduleupdate(::MultilevelBranchingInput) = :branching
getmoduleupdate(::MultilevelBranchingMoranInput) = :moran
getmoduleupdate(::MultilevelMoranInput) = :moran

"""
    update_population_neutral!(population, birthrate, deathrate, moranrate, branchrate, modulesize, 
        branchinitsize, nextID, t, tmax, rng)

Perform a single Gillespie step to advance the simulation of a multilevel population, and
increase the time t → t + Δt.

Calculates the transition rates for a Moran, birth, death or module branching event; then 
selects from these transitions with probability proportional to rate. Time is increased by
`Δt ~ Exp(sum(transitionrates))`. If new time exceeds `tmax` do not perform any transition.
"""

function update_population_neutral!(population, transitionrates, branchrate, modulesize, 
    branchinitsize, modulebranching, t, nextID, nextmoduleID, μ, mutationdist, tmax, 
    maxmodules, moranincludeself, rng; moduleupdate=:branching, timefunc=exptime)

    t += timefunc(rng, sum(transitionrates))
    #only update the population if t < tmax
    if t < tmax
        #choose transition type: 1=moran, 2=asymmetric, 3=birth, 4=death, 5=branch
        transitionid = sample(
            rng, 
            1:5, 
            ProbabilityWeights(transitionrates ./ sum(transitionrates))
        )
        population, nextID, nextmoduleID = transition!(
            population, 
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
        if transitionid in 3:5
            update_neutral_transitionrates!(transitionrates, population, branchrate, modulesize)
        end
    end
    return population, transitionrates, t, nextID, nextmoduleID
end

"""
    transition!(population, transitionid, modulesize, branchinitsize, 
    modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
    moranincludeself, rng; moduleupdate=:branching)
        μ, rng)
    
Perform a single transition step on `population`, determined by `transitionid`.
"""
function transition!(population, transitionid, modulesize, branchinitsize, 
    modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
    moranincludeself, rng; moduleupdate=:branching)
    
    if transitionid == 1
        _, nextID = moranupdate!(population, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself)
    elseif transitionid == 2
        _, nextID = asymmetricupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
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

"""
    moranupdate!(population, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself=false)

Selects a homeostatic module uniformly at random to undergo a single Moran update. From 
that module one cell divides and one cell dies.
"""
function moranupdate!(population, modulesize, t, nextID, μ, mutationdist, rng; moranincludeself=true)
    homeostaticmoduleid, parentcellid, deadcellid = 
        choose_homeostaticmodule_cells(population, rng; moranincludeself, maxmodulesize=modulesize)
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
        _, _, nextID = celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, μ, mutationdist, rng)
    celldeath!(homeostaticmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
    updatetime!(homeostaticmodule, t)
    return population, nextID
end

"""
    asymmetricupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a homeostatic module uniformly at random to undergo a single asymmetric update. From 
that module one cell divides, producing a single offspring.
"""
function asymmetricupdate!(population, modulesize, t, nextID, μ, mutationdist, rng)
    homeostaticmoduleid, parentcellid =
        choose_homeostaticmodule_cells(population, rng; twocells=false, maxmodulesize=modulesize)
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
    _, _, nextID = celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, μ, mutationdist, rng; nchildcells=1)
    updatetime!(homeostaticmodule, t)
    return population, nextID
end


"""
    birthupdate!(population, maxmodulesize, t, nextID, μ, mutationdist, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate!(population, maxmodulesize, t, nextID, μ, mutationdist, rng)
    growingmoduleid, parentcellid = choose_growingmodule_cell(population, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    _, _, nextID = celldivision!(growingmodule, population.subclones, parentcellid, t, nextID, μ, mutationdist, rng)
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
    parentmodule = population.homeostatic_modules[parentmoduleid]
    parentmodule, newmodule, nextID = 
        newmoduleformation!(parentmodule, population.subclones, nextmoduleID, branchinitsize, t, rng; 
            modulebranching, nextID, μ, mutationdist)
    push!(population.growing_modules, newmodule)
    if modulebranching == :split
        move_module_to_growing!(population, parentmoduleid)
    end
    return population, nextmoduleID + 1, nextID
end

function modulemoranupdate!(population, nextmoduleID, branchinitsize, t, rng; 
    modulebranching=:split, nextID=nothing, μ=nothing, mutationdist=nothing)

    population, = modulebranchingupdate!(population, nextmoduleID, branchinitsize, t, rng; 
        modulebranching, nextID, μ, mutationdist)
    deadmoduleid = choose_any_module(population, rng)
    moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
    return population, nextmoduleID + 1, nextID
end

function newmoduleformation!(parentmodule, subclones, nextmoduleID, branchinitsize::Int, t, rng; 
        modulebranching=:split, nextID=nothing, μ=nothing, mutationdist=nothing)
    if modulebranching == :withreplacement
        return sample_new_module_with_replacement!(parentmodule, subclones, nextmoduleID, 
            branchinitsize, t, nextID, μ, mutationdist, rng)
    elseif modulebranching == :withoutreplacement
        return sample_new_module_without_replacement!(parentmodule, subclones, nextmoduleID, 
            branchinitsize, t, nextID, μ, mutationdist, rng)
    elseif modulebranching == :withreplacement_nomutations
                return sample_new_module_with_replacement!(parentmodule, subclones, nextmoduleID, 
                    branchinitsize, t, nextID, 0, :fixed, rng)
    elseif modulebranching == :withoutreplacement_nomutations
        return sample_new_module_without_replacement!(parentmodule, subclones, nextmoduleID, 
            branchinitsize, t, nextID, 0, :fixed, rng)
    elseif modulebranching == :split
        cellmodule, newcellmodule =
            sample_new_module_split!(parentmodule, nextmoduleID, 
                branchinitsize, t, rng)
        return cellmodule, newcellmodule, nextID
    else error("$modulebranching is not a valid module branching option")
    end
end

"""
    sample_new_module_with_replacement!(cellmodule, subclones, newmoduleid, branchinitsize, branchtime, 
        rng::AbstractRNG)

Sample cells with replacement, uniformly at random, from `cellmodule`. Each of these cells 
divides, with one offspring remaining in `cellmodule`, the other becoming part of the new 
module. 

The number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module_with_replacement!(cellmodule, subclones, nextmoduleID, branchinitsize, 
    branchtime, nextID, μ, mutationdist, rng::AbstractRNG)

    ncells = length(cellmodule.cells)
    newcells = eltype(cellmodule.cells)[]
    for i in 1:branchinitsize
        randcellid = rand(rng, 1:ncells)
        cellmodule, subclones, nextID = 
            celldivision!(cellmodule, subclones, randcellid, branchtime, nextID, μ, mutationdist, rng)
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
    sample_new_module_without_replacement!(cellmodule::T, subclones, nextmoduleID, branchinitsize, 
    branchtime, nextID, μ, mutationdist, rng::AbstractRNG) where T<: AbstractModule

Sample cells without replacement, uniformly at random, from `cellmodule`. Each of these 
cells divides, with one offspring remaining in `cellmodule`, the other becoming part of the 
new module. 

The number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module_without_replacement!(cellmodule, subclones, nextmoduleID, branchinitsize, 
    branchtime, nextID, μ, mutationdist, rng::AbstractRNG)

    sampleids = sample(rng, 1:length(cellmodule.cells), branchinitsize, replace=false)
    newcells = eltype(cellmodule.cells)[]
    for cellid in sampleids
        cellmodule, subclones, nextID = 
            celldivision!(cellmodule, subclones, cellid, branchtime, nextID, μ, mutationdist, rng)
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
    sample_new_module_split!(cellmodule::T, nextmoduleID, branchinitsize, branchtime, 
    rng::AbstractRNG) where T<: AbstractModule

Sample cells without replacement, uniformly at random, from `cellmodule` to form the new 
module. Remaining cells form the parent module.

The number of cells is given by `branchinitsize` and the newmodule is given initial time
`branchtime`.

"""
function sample_new_module_split!(cellmodule, nextmoduleID, branchinitsize, branchtime, 
    rng::AbstractRNG)

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
    choose_any_module(population, rng::AbstractRNG)
Select a module uniformly at random from all modules. Returns `moduleid` and the vector it
is in (either `growing_modules` or `homeostatic_modules`.)
"""
function choose_any_module(population, rng::AbstractRNG)
    moduleid = rand(rng, 1:length(population))
    return moduleid
end

"""
    choose_homeostaticmodule_cells(population, maxmodulesize, rng::AbstractRNG)
Select a homeostatic module and the ids of two cells from the module, uniformly at random.
"""
function choose_homeostaticmodule_cells(population, rng::AbstractRNG; twocells=true, 
    moranincludeself=true, maxmodulesize=length(first(population.homeostatic_modules)))

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
    get_neutral_transitionrates(population, birthrate, deathrate, moranrate, branchrate, modulesize)
    Compute the rates for moran update, asymmetric update, birth update, death update and
        module branching.
"""
function get_neutral_transitionrates(population, branchrate, modulesize)
    rates = zeros(Float64, 5)
    return update_neutral_transitionrates!(rates, population, branchrate, modulesize)
end

function update_neutral_transitionrates!(rates, population, branchrate, modulesize)
    birthrate, deathrate, moranrate, asymmetricrate = getwildtyperates(population)
    number_homeostatic_modules = length(population.homeostatic_modules)
    number_cells_in_growing_modules = sum(length.(population.growing_modules))
    rates[1] = number_homeostatic_modules * moranrate * modulesize
    rates[2] = number_homeostatic_modules * asymmetricrate * modulesize
    rates[3] = number_cells_in_growing_modules * birthrate
    rates[4] = number_cells_in_growing_modules * deathrate
    rates[5] = number_homeostatic_modules * branchrate
    return rates
end

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

function removemodule!(population, dyingmoduleid; moduletype=:all)
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

function move_module_to_homeostasis!(population, cellmoduleid::Integer)
    cellmodule = popat!(population.growing_modules, cellmoduleid)
    push!(population.homeostatic_modules, cellmodule)
    return population
end

function move_module_to_growing!(population, cellmoduleid::Integer)
    cellmodule = popat!(population.homeostatic_modules, cellmoduleid)
    push!(population.growing_modules, cellmodule)
    return population
end

killallcells!(population::CellVector, args...) = nothing