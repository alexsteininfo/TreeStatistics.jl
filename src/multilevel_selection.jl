
function simulate!(population, selection, tmax, maxmodules, branchrate, 
    modulesize, branchinitsize, modulebranching, μ, mutationdist, moranincludeself, nextID, 
    nextmoduleID, rng; moduleupdate=:branching, t0=nothing)

    t = isnothing(t0) ? age(population) : t0
    nsubclones = getmaxsubclones(selection)
    nsubclonescurrent = length(population.subclones)
    transitionrates = get_selection_transitionrates(population, branchrate, nsubclones)

    while t < tmax && (moduleupdate==:moran || length(population) < maxmodules)

        population, transitionrates, nsubclonescurrent, t, nextID, nextmoduleID = 
            update_population_selection!(population, transitionrates, nsubclonescurrent, 
                nsubclones, selection, branchrate, modulesize, 
                branchinitsize, modulebranching, t, nextID, nextmoduleID, μ, mutationdist, 
                tmax, maxmodules, moranincludeself, rng; moduleupdate)

        #returns empty list of modules if population dies out
        if length(population) == 0
            break
        end
    end
    return population, nextID, nextmoduleID
end


function update_population_selection!(population, transitionrates, nsubclonescurrent, 
    nsubclones, selection, branchrate, modulesize, branchinitsize, modulebranching, 
    t, nextID, nextmoduleID, μ, mutationdist, tmax, maxmodules, moranincludeself, rng; 
    moduleupdate=:branching)

    t += exptime(rng, sum(transitionrates))
    #only update the population if t < tmax
    
    if t < tmax
        #choose transition type: 
        #1 0*nsubclones + (1:nsubclones) = moran_1:moran_nsubclones,
        #2 1*nsubclones + (1:nsubclones) = asymmetric_1:asymmetric_nsubclones,
        #3 2*nsubclones + (1:nsubclones) = birth_1:birth_nsubclones,
        #4 3*nsubclones + (1:nsubclones) = death_1:death_nsubclones,
        #5 4*nsubclones + 1 = moduleformation
        transitionid = sample(
            rng, 
            1:(5*nsubclones+1), 
            ProbabilityWeights(transitionrates ./ sum(transitionrates))
        )
        transition_type = (transitionid - 1) ÷ nsubclones + 1 # in 1:5 as above
        transition_subcloneid = (transitionid - 1) % nsubclones + 1 #id of subclone

        population, transitionrates, nsubclonescurrent, nextID, nextmoduleID = transition_selection!(
            population, 
            transitionrates,
            transition_type, 
            transition_subcloneid,
            nsubclonescurrent,
            nsubclones,
            selection,
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
            branchrate,
            rng;
            moduleupdate
        )
    end
    return population, transitionrates, nsubclonescurrent, t, nextID, nextmoduleID
end

"""
    transition_selection!(population, transitionrates, transition_type, transition_subcloneid, nsubclonescurrent,
    nsubclones, mutant_selection, mutant_time, modulesize, branchinitsize, 
    modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
    moranincludeself, branchrate, rng; moduleupdate=:branching)
    
Perform a single transition step on `population`, determined by `transition_type` and
`transition_subclone`.
"""
function transition_selection!(population, transitionrates, transition_type, 
    transition_subcloneid, nsubclonescurrent, nsubclones, selection, modulesize, 
    branchinitsize, modulebranching, t, nextID, nextmoduleID, μ, mutationdist, maxmodules, 
    moranincludeself, branchrate, rng; moduleupdate=:branching)
    
    if transition_type == 1
        population, transitionrates, nsubclonescurrent, nextID = moranupdate_selection!(
            population, 
            transitionrates,
            transition_subcloneid, 
            nsubclonescurrent,
            nsubclones, 
            selection,
            modulesize, 
            t, 
            nextID, 
            μ, 
            mutationdist, 
            rng; 
            moranincludeself
        )
    elseif transition_type == 2
        population, transitionrates, nsubclonescurrent, nextID = asymmetricupdate_selection!(
            population, 
            transitionrates,
            transition_subcloneid, 
            nsubclonescurrent,
            nsubclones, 
            selection,
            t, 
            nextID, 
            μ, 
            mutationdist, 
            rng
        )
    elseif transition_type == 3
        population, transitionrates, nsubclonescurrent, nextID = birthupdate_selection!(
            population, 
            transitionrates,
            transition_subcloneid, 
            nsubclonescurrent,
            nsubclones, 
            selection,
            modulesize, 
            t, 
            nextID, 
            μ, 
            mutationdist, 
            branchrate,
            rng
        )
    elseif transition_type == 4
        population, transitionrates = deathupdate_selection!(
            population, 
            transitionrates, 
            transition_subcloneid, 
            nsubclones, 
            t, 
            μ, 
            mutationdist, 
            rng
        )
    elseif transition_type == 5
        if moduleupdate == :branching || length(population) < maxmodules
            population, transitionrates, nextmoduleID, nextID = modulebranchingupdate_selection!(
                population, 
                transitionrates,
                nextmoduleID, 
                branchinitsize, 
                nsubclones,
                modulesize,
                t, 
                branchrate,
                rng; 
                modulebranching, 
                nextID, 
                μ,
                mutationdist
            )
        elseif moduleupdate == :moran
            population, transitionrates, nextmoduleID, nextID = modulemoranupdate_selection!(
                population, 
                transitionrates,
                nextmoduleID, 
                branchinitsize, 
                nsubclones,
                modulesize,
                t, 
                branchrate,
                rng; 
                modulebranching, 
                nextID, 
                μ,
                mutationdist
            )
        end
    end
    return population, transitionrates, nsubclonescurrent, nextID, nextmoduleID
end

"""
    moranupdate_selection!(population, transitionrates, transition_subcloneid, nsubclonescurrent,
    nsubclones, mutant_selection, mutant_time, modulesize, t, nextID, μ, mutationdist, rng; 
    moranincludeself=true)

Selects a homeostatic module uniformly at random to undergo a single Moran update. From 
that module one cell divides and one cell dies.
"""
function moranupdate_selection!(population, transitionrates, transition_subcloneid, nsubclonescurrent,
    nsubclones, selection::AbstractSelection, modulesize, t, nextID, μ, mutationdist, rng; 
    moranincludeself=true)
    #choose cells to die and divide
    moranrate = population.subclones[transition_subcloneid].moranrate
    N_transition_cells = round(Int64, (transitionrates[transition_subcloneid]) / moranrate) 
    homeostaticmoduleid, parentcellid = 
        choose_module_cell(population.homeostatic_modules, transition_subcloneid, 
            N_transition_cells, rng)

    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
    deadcellid = 
        choose_moran_deadcell(modulesize, parentcellid, moranincludeself, rng)

    #implement cell division
    homesotaticmodule, subclones, nextID = 
        celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, 
            μ, mutationdist, rng)
    parentcell = homeostaticmodule[parentcellid]
    deadcell = homeostaticmodule[deadcellid]
    celldeath!(homeostaticmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)

    #update transition rates if necessary
    transitionrates = 
        update_selection_transitionrates_after_moran!(
            transitionrates, population, getclonetype(parentcell), 
                getclonetype(deadcell), nsubclones)   

    #check if ready for a new mutant subclone
    if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
        population, homeostaticmodule, nsubclonescurrent, transitionrates = 
            update_cellmutation!(population, homeostaticmodule, nsubclonescurrent, nsubclones,
                transitionrates, parentcell, selection, t, :homeostatic, rng)
    end
    updatetime!(homeostaticmodule, t)
    return population, transitionrates, nsubclonescurrent, nextID
end

"""
    asymmetricupdate_selection!(population, modulesize, t, nextID, μ, mutationdist, rng)

Selects a homeostatic module uniformly at random to undergo a single asymmetric update. From 
that module one cell divides, producing a single offspring.
"""
function asymmetricupdate_selection!(population, transitionrates, transition_subcloneid, nsubclonescurrent,
    nsubclones, selection::AbstractSelection, t, nextID, μ, mutationdist, rng)
    #choose cell to divide
    asymmetricrate = population.subclones[transition_subcloneid].asymmetricrate
    N_transition_cells = round(Int64, (transitionrates[nsubclones + transition_subcloneid]) / asymmetricrate) 
    homeostaticmoduleid, parentcellid = 
        choose_module_cell(population.homeostatic_modules, transition_subcloneid, 
            N_transition_cells, rng)
    homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
    parentcell = homeostaticmodule[parentcellid]

    homeostaticmodule, subclones, nextID = 
        celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, μ, 
            mutationdist, rng; nchildcells=1)
    #check if ready for a new mutant subclone
    if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
        population, homeostaticmodule, nsubclonescurrent, transitionrates = 
            update_cellmutation!(population, homeostaticmodule, nsubclonescurrent, nsubclones,
                transitionrates, parentcell, selection, t, :homeostatic, rng)
    end
    updatetime!(homeostaticmodule, t)
    return population, transitionrates, nsubclonescurrent, nextID
end


"""
    birthupdate_selection!(population, transitionrates, transition_subcloneid, nsubclonescurrent,
    nsubclones, mutant_selection, mutant_time, maxmodulesize, t, nextID, μ, mutationdist, branchrate, rng; 
    moranincludeself=true)

Selects a cell uniformly at random from all cells in non-homeostatic modules to divide.
"""
function birthupdate_selection!(population, transitionrates, transition_subcloneid, nsubclonescurrent,
    nsubclones, selection::AbstractSelection, maxmodulesize, t, nextID, μ, mutationdist, branchrate, rng)

    #choose cell to divide
    birthrate = population.subclones[transition_subcloneid].birthrate
    N_transition_cells = round(Int64, (transitionrates[2*nsubclones + transition_subcloneid]) / birthrate) 
    growingmoduleid, parentcellid = 
        choose_module_cell(population.growing_modules, transition_subcloneid, 
            N_transition_cells, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    parentcell = growingmodule[parentcellid]

    transitiontohomeostasis = 
        if length(growingmodule) + 1 == maxmodulesize
            move_module_to_homeostasis!(population, growingmoduleid)
            true
        else
            false
        end 

    growingmodule, subclones, nextID = celldivision!(growingmodule, population.subclones, 
        parentcellid, t, nextID, μ, mutationdist, rng)

    #update transition rates if necessary
    transitionrates = 
        update_selection_transitionrates_after_birth!(
            transitionrates, population, getclonetype(parentcell), 
                nsubclones, transitiontohomeostasis, branchrate, growingmodule)  
    updatetime!(growingmodule, t)
    
    #check if ready for a new mutant subclone
    if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
        moduletype = transitiontohomeostasis ? :homeostatic : :growing
        population, growingmodule, nsubclonescurrent, transitionrates = 
            update_cellmutation!(population, growingmodule, nsubclonescurrent, nsubclones,
                transitionrates, parentcell, selection, t, moduletype, rng)
    end
    return population, transitionrates, nsubclonescurrent, nextID
end

"""
    deathupdate_selection!(population, modulesize, t, rng)

Selects a cell uniformly at random from all cells in non-homeostatic modules to die. If 
cell death results in an empty module, remove that module from the population.
"""

function deathupdate_selection!(population, transitionrates, transition_subcloneid, 
    nsubclones, t, μ, mutationdist, rng)

    deathrate = population.subclones[transition_subcloneid].deathrate
    N_transition_cells = round(Int64, (transitionrates[3*nsubclones + transition_subcloneid]) / deathrate) 
    growingmoduleid, deadcellid = 
        choose_module_cell(population.growing_modules, transition_subcloneid, 
        N_transition_cells, rng)
    growingmodule = population.growing_modules[growingmoduleid]
    deadcell = growingmodule[deadcellid]
    moduleextinct = length(growingmodule) == 1
    celldeath!(growingmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
    transitionrates = 
        update_selection_transitionrates_after_death!(transitionrates, population, 
            getclonetype(deadcell), nsubclones)   
    updatetime!(growingmodule, t)
    if moduleextinct
        moduledeath!(population, growingmoduleid; moduletype=:growing)
    end
    return population, transitionrates
end

function modulebranchingupdate_selection!(population, transitionrates, nextmoduleID, 
    branchinitsize, nsubclones, maxmodulesize, t, branchrate, rng; modulebranching=:split, nextID=nothing, 
    μ=nothing, mutationdist=nothing)
    
    parentmoduleid = choose_homeostaticmodule(population, rng)
    parentmodule = population.homeostatic_modules[parentmoduleid]
    parentmodule, newmodule, nextID = 
        newmoduleformation!(parentmodule, population.subclones, nextmoduleID, branchinitsize, t, rng; 
            modulebranching, nextID, μ, mutationdist)
    push!(population.growing_modules, newmodule)
    if modulebranching == :split
        move_module_to_growing!(population, parentmoduleid)
    end
    transitionrates = update_selection_transitionrates_after_newmodule!(
        transitionrates, population, parentmodule, newmodule, nsubclones, branchrate, maxmodulesize)
    return population, transitionrates, nextmoduleID + 1, nextID
end

function modulemoranupdate_selection!(population, transitionrates, nextmoduleID, 
    branchinitsize, nsubclones, maxmodulesize, t, branchrate, rng; modulebranching=:split, nextID=nothing, 
    μ=nothing, mutationdist=nothing)

    population, = modulebranchingupdate_selection!(population,transitionrates, nextmoduleID, 
        branchinitsize, nsubclones, maxmodulesize, t, branchrate, rng; modulebranching, nextID, μ, 
        mutationdist)
    deadmoduleid = choose_any_module(population, rng)
    deadmodule = population[deadmoduleid]
    moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
    transitionrates = update_selection_transitionrates_after_moduledeath!(
        transitionrates, population, deadmodule, nsubclones, branchrate, maxmodulesize)
    return population, transitionrates, nextmoduleID + 1, nextID
end

function update_cellmutation!(population, mutatingmodule, nsubclonescurrent, nsubclones,
    transitionrates, parentcell, selection::AbstractSelection, t, moduletype, rng)

    old_clonetype = getclonetype(parentcell)
    new_clonetype = nsubclonescurrent + 1
    newmutant_selection = getselectioncoefficient(selection, nsubclonescurrent, rng)
    cellmutation!(mutatingmodule, population.subclones, newmutant_selection, parentcell, t)
    transitionrates = update_selection_transitionrates_after_newsubclone!(
            transitionrates, population, new_clonetype, old_clonetype, nsubclones, moduletype)   
    nsubclonescurrent += 1
    return population, mutatingmodule, nsubclonescurrent, transitionrates
end

"""
    get_selection_transitionrates(population, branchrate, modulesize, nsubclones)
Compute the rates for moran update, asymmetric update, birth update and death update for 
    each subclone and module formation rate and return as a Vector in the form
     `[moran_1..., moran_nsubclones, asymmetric_1..., birth_1..., death1..., death_nsubclones,
     moduleformation]`
"""
function get_selection_transitionrates(population, branchrate, nsubclones)
    rates = zeros(Float64, 4 * nsubclones + 1)
    return update_selection_transitionrates!(rates, population, branchrate, nsubclones)
end

function update_selection_transitionrates!(rates, population, branchrate, nsubclones)
    nsubclonescurrent = length(population.subclones)
    ncells_homeostatic_by_subclone = 
        number_cells_by_subclone(population.homeostatic_modules, nsubclonescurrent)
    ncells_growing_by_subclone = 
        number_cells_by_subclone(population.growing_modules, nsubclonescurrent)

    for (i, (ncells_homeostatic, ncells_growing, subclone)) in enumerate(zip(
        ncells_homeostatic_by_subclone, 
        ncells_growing_by_subclone, 
        population.subclones
    ))
        rates[i + 0 * nsubclones] = ncells_homeostatic * subclone.moranrate
        rates[i + 1 * nsubclones] = ncells_homeostatic * subclone.asymmetricrate
        rates[i + 2 * nsubclones] = ncells_growing * subclone.birthrate
        rates[i + 3 * nsubclones] = ncells_growing * subclone.deathrate
    end
    rates[end] = length(population.homeostatic_modules) * branchrate
    return rates
end

function update_selection_transitionrates_after_moran!(rates, population, divided_subcloneid, dead_subcloneid, nsubclones)
    if divided_subcloneid == dead_subcloneid
        return rates
    else
        #update moranrates
        rates[divided_subcloneid] += population.subclones[divided_subcloneid].moranrate
        rates[dead_subcloneid] -= population.subclones[dead_subcloneid].moranrate
        #update asymmetricrates
        rates[nsubclones + divided_subcloneid] += population.subclones[divided_subcloneid].asymmetricrate
        rates[nsubclones + dead_subcloneid] -= population.subclones[dead_subcloneid].asymmetricrate
    end
    return rates
end

function update_selection_transitionrates_after_newsubclone!(rates, population, new_subcloneid, old_subcloneid, nsubclones, moduletype)
    if moduletype == :homeostatic
        #update moranrates
        rates[new_subcloneid] += population.subclones[new_subcloneid].moranrate
        rates[old_subcloneid] -= population.subclones[old_subcloneid].moranrate
        #update asymmetricrates
        rates[nsubclones + new_subcloneid] += population.subclones[new_subcloneid].asymmetricrate
        rates[nsubclones + old_subcloneid] -= population.subclones[old_subcloneid].asymmetricrate
    elseif moduletype == :growing
        #update birthrates
        rates[2 * nsubclones + new_subcloneid] += population.subclones[new_subcloneid].birthrate
        rates[2 * nsubclones + old_subcloneid] -= population.subclones[old_subcloneid].birthrate
        #update deathrates
        rates[3 * nsubclones + new_subcloneid] += population.subclones[new_subcloneid].deathrate
        rates[3 * nsubclones + old_subcloneid] -= population.subclones[old_subcloneid].deathrate
    else error("$moduletype not allowed `moduletype`")
    end
    return rates
end

function update_selection_transitionrates_after_birth!(rates, population, divided_subcloneid, nsubclones, transitiontohomeostasis, branchrate, transformedmodule)
    if transitiontohomeostasis
        #move original cell rates from growing to homeostatic processes
        for cell in transformedmodule.cells
            cell_subclone = population.subclones[getclonetype(cell)]
            rates[getclonetype(cell)] += cell_subclone.moranrate
            rates[nsubclones + getclonetype(cell)] += cell_subclone.asymmetricrate
            rates[2 * nsubclones + getclonetype(cell)] -= cell_subclone.birthrate
            rates[3 * nsubclones + getclonetype(cell)] -= cell_subclone.deathrate
        end
        #correct above for new cell (was not included in birth/death rates)
        rates[2 * nsubclones + divided_subcloneid] += population.subclones[divided_subcloneid].birthrate
        rates[3 * nsubclones + divided_subcloneid] += population.subclones[divided_subcloneid].deathrate
        #additional homeostatic module so update branchrate
        rates[end] += branchrate
    else
        rates[2 * nsubclones + divided_subcloneid] += population.subclones[divided_subcloneid].birthrate
        rates[3 * nsubclones + divided_subcloneid] += population.subclones[divided_subcloneid].deathrate
    end
    return rates
end

function update_selection_transitionrates_after_death!(rates, population, dead_subcloneid, nsubclones)

    #update birthrates
    rates[2 * nsubclones + dead_subcloneid] -= population.subclones[dead_subcloneid].birthrate
    #update deathrates
    rates[3 * nsubclones + dead_subcloneid] -= population.subclones[dead_subcloneid].deathrate
    return rates
end

function update_selection_transitionrates_after_newmodule!(rates, population, 
    parentmodule, newmodule, nsubclones, branchrate, maxmodulesize)
    
    #add to birth/death rates for new module cells
    #if length(parentmodule) < maxmodulesize these cells have been removed from a
    #homeostatic module so also update moran and asymmetric rates
    for cell in newmodule.cells
        cell_subclone = population.subclones[getclonetype(cell)]
        rates[2 * nsubclones + getclonetype(cell)] += cell_subclone.birthrate
        rates[3 * nsubclones + getclonetype(cell)] += cell_subclone.deathrate
        if length(parentmodule) < maxmodulesize
            rates[getclonetype(cell)] -= cell_subclone.moranrate
            rates[nsubclones + getclonetype(cell)] -= cell_subclone.asymmetricrate
        end
    end

    if length(parentmodule) < maxmodulesize
        #parent module gone from homeostatic to growing
        for cell in parentmodule.cells
            cell_subclone = population.subclones[getclonetype(cell)]
            rates[getclonetype(cell)] -= cell_subclone.moranrate
            rates[nsubclones + getclonetype(cell)] -= cell_subclone.asymmetricrate
            rates[2 * nsubclones + getclonetype(cell)] += cell_subclone.birthrate
            rates[3 * nsubclones + getclonetype(cell)] += cell_subclone.deathrate
        end
        rates[end] -= branchrate
    end

    return rates
end

function update_selection_transitionrates_after_moduledeath!(rates, population, 
    dyingmodule, nsubclones, branchrate, maxmodulesize)

    if length(dyingmodule) < maxmodulesize
        for cell in dyingmodule.cells
            cell_subclone = population.subclones[getclonetype(cell)]
            rates[2 * nsubclones + getclonetype(cell)] -= cell_subclone.birthrate
            rates[3 * nsubclones + getclonetype(cell)] -= cell_subclone.deathrate
        end
    elseif length(dyingmodule) == maxmodulesize
        for cell in dyingmodule.cells
            cell_subclone = population.subclones[getclonetype(cell)]
            rates[getclonetype(cell)] -= cell_subclone.moranrate
            rates[nsubclones + getclonetype(cell)] -= cell_subclone.asymmetricrate
        end
        rates[end] -= branchrate
    else 
        error("$modulesize not valid modulesize")
    end
    return rates
end
"""
    choose_module_cell(modules, subclone, rng::AbstractRNG)
Select a cell uniformly at random from all cells in `modules `in the
given `subclone`, and return the module and cell id.
"""
# function choose_module_cell(modules, subclone, rng::AbstractRNG)
#     subclone_module_cell_ids = get_module_cells_given_subclones(modules, subclone)
#     moduleid, cellid = rand(
#         rng, 
#         subclone_module_cell_ids
#     )
#     return moduleid, cellid
# end
#TODO implement this method properly for all update types. Should I instead keep track of a 
#list of cells in each compartment??
function choose_module_cell(modules, subclone, N_transition_cells, rng::AbstractRNG)
    chosen_cell = rand(rng, 1:N_transition_cells)
    return get_moduleid_cellid(modules, subclone, chosen_cell)
end

function get_moduleid_cellid(modules, subclone, chosen_cell)
    i = 1
    for (moduleid, mod) in enumerate(modules)
        for (cellid, cell) in enumerate(mod.cells)
            if cell.clonetype == subclone
                chosen_cell == i && return moduleid, cellid
                i +=1
            end
        end
    end
end

function allclonetypes(modules)
    return [getclonetype(cell) for mod in modules for cell in mod.cells]
end

function get_module_cells_given_subclones(modules, subclone)
    return [(moduleid, cellid) 
        for (moduleid, mod) in enumerate(modules)
            for (cellid, cell) in enumerate(mod.cells)
                if cell.clonetype == subclone
    ]
end
