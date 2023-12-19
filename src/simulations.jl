function simulate!(population::SinglelevelPopulation, input::BranchingMoranInput, 
    selection::AbstractSelection, nextID::Int64, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, t0=nothing, tmax=nothing)

    population, nextID = branchingprocess!(
        population,
        selection,
        input.Nmax, 
        input.μ, 
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng;
        timefunc,
        t0
    )
    
    population, nextID = moranprocess!(
        population,
        selection,
        input.μ, 
        input.mutationdist, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng;
        timefunc,
        t0,
        moranincludeself=input.moranincludeself
    )

    return population, nextID
end

function simulate!(population::SinglelevelPopulation, input::BranchingInput, 
    selection::AbstractSelection, nextID::Int64, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, t0=nothing, tmax=nothing)

    population, nextID = branchingprocess!(
        population,
        selection,
        input.Nmax, 
        input.μ, 
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng; 
        timefunc,
        t0
    )
    
    return population, nextID
end

function simulate!(population::SinglelevelPopulation, input::MoranInput, 
    selection::AbstractSelection, nextID::Int64, rng::AbstractRNG=Random.GLOBAL_RNG; 
    timefunc=exptime, t0=nothing, tmax=nothing)

    population, nextID = moranprocess!(
        population,
        selection,
        input.μ, 
        input.mutationdist, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng; 
        timefunc,
        t0,
        moranincludeself=input.moranincludeself
    )

    return population, nextID
end

"""
    branchingprocess!(population::SinglelevelPopulation{CellModule}, 
    selection::AbstractSelection, Nmax, μ, mutationdist, tmax, nextID, rng::AbstractRNG;
    timefunc=exptime, t0=nothing)

Run branching process simulation, starting in state defined by cellmodule.

Simulate a stochastic branching process, starting with a single cell, with with birth rate 
`b`, death rate `d` until population reaches size `Nmax`.

Cells accumulate neutral mutations at division with rate `μ`.
"""
function branchingprocess!(population::SinglelevelPopulation, selection::AbstractSelection, 
    Nmax, μ, mutationdist, tmax, nextID, rng::AbstractRNG; timefunc=exptime, t0=nothing)

    t = !isnothing(t0) ? t0 : age(population.singlemodule)
    N = length(population.singlemodule)

    nsubclones = getmaxsubclones(selection)
    nsubclonescurrent = length(population.subclones)
    birthrates = getbirthrates(population.subclones)
    deathrates = getdeathrates(population.subclones)

    #Rmax starts with birthrate + deathrate and changes once a fitter mutant is introduced
    Rmax = maximum(birthrates) + maximum(deathrates)

    while N < Nmax && N > 0

        #calc next event time and break if it exceeds tmax 
        Δt =  1 / (Rmax * N) .* timefunc(rng)
        t + Δt <= tmax || break # end simulation if time exceeds maximum
        t += Δt

        population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones = 
            branchingupdate!(population, selection, birthrates, deathrates, Rmax, N, nextID, 
                nsubclonescurrent, nsubclones, t, μ, 
                mutationdist, rng)
    end
    return population, nextID
end

function branchingupdate!(population::SinglelevelPopulation, selection, birthrates, 
    deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones, t, μ, mutationdist, rng)

    randcellid = rand(rng, 1:N) #pick a random cell
    randcell = population.singlemodule.cells[randcellid]
    r = rand(rng, Uniform(0, Rmax))
    #get birth and death rates for randcell
    cellsubclone = getclonetype(randcell)
    br = birthrates[cellsubclone]
    dr = deathrates[cellsubclone]

    if r < br 
        #cell divides
        cellmodule, subclones, nextID = celldivision!(population.singlemodule, population.subclones, 
            randcellid, t, nextID, μ, mutationdist, rng)
        N += 1
        #check if ready for a new mutant subclone
        if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
            newmutant_selectioncoeff = getselectioncoefficient(selection, nsubclonescurrent, rng)
            cellmutation!(cellmodule, population.subclones, newmutant_selectioncoeff, randcell, t)
            nsubclonescurrent += 1
            #update rates and Rmax
            birthrates = getbirthrates(population.subclones)
            deathrates = getdeathrates(population.subclones)
            Rmax = maximum(birthrates) + maximum(deathrates)
        end
    elseif r < br + dr
        #cell dies
        cellmodule = celldeath!(population.singlemodule, population.subclones, randcellid, t)
        N -= 1
        #return empty population if all cells have died
    end
    updatetime!(population.singlemodule, t)
    return population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones
end

"""
    moranprocess!(population::SinglelevelPopulation, tmax, μ, mutationdist, 
    rng::AbstractRNG)

Run Moran process simulation, starting in state defined by cellmodule.
"""
function moranprocess!(population::SinglelevelPopulation, selection, μ, mutationdist, tmax, 
    nextID, rng::AbstractRNG; timefunc=exptime, t0=nothing, moranincludeself=true)

    t = !isnothing(t0) ? t0 : age(population.singlemodule)
    N = length(population.singlemodule)
    nsubclones = getmaxsubclones(selection)
    nsubclonescurrent = length(population.subclones)
    moranrates = getmoranrates(population.subclones)
    Rmax = maximum(moranrates)
    while true

        #calc next event time and break if it exceeds tmax
        Δt =  1 / (Rmax * N) .* timefunc(rng)
        t = t + Δt
        if t > tmax
            break
        end
        population, moranrates, Rmax, N, nextID, nsubclonescurrent = 
            moranupdate!(population, selection, moranrates, Rmax, N, nextID, 
                nsubclonescurrent, nsubclones, t, μ, mutationdist, moranincludeself, rng)

    end
    return population, nextID
end

function moranupdate!(population::SinglelevelPopulation, selection, moranrates, Rmax, N, nextID, 
    nsubclonescurrent, nsubclones, t, μ, mutationdist, moranincludeself, rng)

    dividecellid = rand(rng, 1:N) #pick a random cell to divide 
    dividecell = population.singlemodule.cells[dividecellid]
    r = rand(rng, Uniform(0, Rmax))
    mr = moranrates[getclonetype(dividecell)]
    if r < mr
        #pick a random cell to die
        deadcell = choose_moran_deadcell(N, dividecellid, moranincludeself, rng)
        #cell divides
        cellmodule, subclones, nextID = celldivision!(population.singlemodule, population.subclones, 
            dividecellid, t, nextID, μ, mutationdist, rng)

        #check if t>=mutant_time for next fit subclone and more subclones are expected
            if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
                newmutant_selectioncoeff = getselectioncoefficient(selection, nsubclonescurrent, rng)
                cellmutation!(cellmodule, population.subclones, newmutant_selectioncoeff, dividecell, t)
            nsubclonescurrent += 1
            #update rates and Rmax
            moranrates = getmoranrates(population.subclones)
            Rmax = maximum(moranrates)
        end

        #cell dies
        cellmodule = celldeath!(cellmodule, population.subclones, deadcell, t)
    end
    updatetime!(population.singlemodule, t)
    return population, moranrates, Rmax, N, nextID, nsubclonescurrent

end

function choose_moran_deadcell(modulesize, dividecellid, moranincludeself, rng)
    if moranincludeself
        deadcellid = rand(rng, 1:modulesize)
        #if dead cell and divide cell are the same kill one of the offspring
        if deadcellid == dividecellid 
            return modulesize + 1 
        else 
            return deadcellid
        end
    else
        #exclude dividecellidx
        return rand(rng, deleteat!(collect(1:modulesize), dividecellid))
    end
end


getbirthrates(subclones) = Float64[subclone.birthrate for subclone in subclones]
getdeathrates(subclones) = Float64[subclone.deathrate for subclone in subclones]
getmoranrates(subclones) = Float64[subclone.moranrate for subclone in subclones]
getasymmetricrates(subclones) = Float64[subclone.asymmetricrate for subclone in subclones]



function getwildtyperates(population::AbstractPopulation)
    return getwildtyperates(population.subclones)
end

function getwildtyperates(subclones::Vector{Subclone})
    return (
        birthrate = subclones[1].birthrate,
        deathrate = subclones[1].deathrate,
        moranrate = subclones[1].moranrate,
        asymmetricrate = subclones[1].asymmetricrate
    )
end

"""
    get_newsubclone_rates(wildtype, selection)
Compute new birth, death, moran and asymmetric rates for a new subclone.

All `wildtype` rates are increased by a factor of `(1 + selection)` (except for death rate 
which is unchanged).
"""
function get_newsubclone_rates(wildtype, selectioncoefficient)
    #TODO check whether this is how Francesco defines (i.e. should selection act on wildtype
    #fitness or parent fitness?)
    return (
        birthrate = wildtype.birthrate * (1 + selectioncoefficient),
        deathrate = wildtype.deathrate,
        moranrate = wildtype.moranrate * (1 + selectioncoefficient),
        asymmetricrate =  wildtype.asymmetricrate * (1 + selectioncoefficient)
    )
end

function addmutations!(cell1::Cell, cell2::Cell, μ, mutID, rng, mutationdist=mutationdist, Δt=Δt)
    if mutationdist == :poissontimedep || mutationdist == :fixedtimedep
        #if mutations are time dependent we add the mutations accumulated by the parent cell
        #to both children at division
        numbermutations = numbernewmutations(rng, mutationdist, μ, Δt=Δt)
        mutID = addnewmutations!(cell1, cell2, numbermutations, mutID)
    else
        numbermutations = numbernewmutations(rng, mutationdist, μ)
        mutID = addnewmutations!(cell1, numbermutations, mutID)
        numbermutations = numbernewmutations(rng, mutationdist, μ)
        mutID = addnewmutations!(cell2, numbermutations, mutID)
    end
    return mutID
end

function addmutations!(cell::Cell, μ, mutID, rng, mutationdist=mutationdist, Δt=Δt)
    if mutationdist == :poissontimedep || mutationdist == :fixedtimedep
        #if mutations are time dependent we add the mutations accumulated by the parent cell
        #child at division
        numbermutations = numbernewmutations(rng, mutationdist, μ, Δt=Δt)
        mutID = addnewmutations!(cell1, numbermutations, mutID)
    else
        numbermutations = numbernewmutations(rng, mutationdist, μ)
        mutID = addnewmutations!(cell, numbermutations, mutID)
    end
    return mutID
end

function numbernewmutations(rng, mutationdist, μ; Δt=nothing)
    if mutationdist == :fixed
        return round(Int64, μ)
    elseif mutationdist == :poisson
        return rand(rng, Poisson(μ))
    elseif mutationdist == :geometric
        return rand(rng, Geometric(1/(1+μ)))
    elseif mutationdist == :poissontimedep
        return rand(rng, Poisson(μ*Δt))
    elseif mutationdist == :fixedtimedep
        return round(Int64, μ*Δt)
    else
        error("$mutationdist is not a valid mutation rule")
    end
end

function addnewmutations!(cell::Cell, numbermutations, mutID)
    #function to add new mutations to cells
    newmutations = mutID:mutID + numbermutations - 1
    append!(cell.mutations, newmutations)
    mutID = mutID + numbermutations
    return mutID
end

function addnewmutations!(cell1::Cell, cell2::Cell, numbermutations, mutID)
    newmutations = mutID:mutID + numbermutations - 1
    append!(cell1.mutations, newmutations)
    append!(cell2.mutations, newmutations)
    mutID = mutID + numbermutations
    return mutID
end

updatetime!(abstractmodule, t) = abstractmodule.t = t

function celldivision!(cellmodule::CellModule, subclones, parentcellid, t, mutID, μ, mutationdist, rng; nchildcells=2)
    
    Δt = t - cellmodule.cells[parentcellid].birthtime
    cellmodule.cells[parentcellid].birthtime = t
    if nchildcells == 2
        push!(cellmodule.cells, copycell(cellmodule.cells[parentcellid])) #add new copy of parent cell to cells
        cellmodule.cells[end].id = cellmodule.cells[end-1].id + 1
        cellmodule.cells[end].parentid = cellmodule.cells[parentcellid].id
        subclones[cellmodule.cells[parentcellid].clonetype].size += 1

    end
    #add new mutations to both new cells
    if μ > 0.0 
        if nchildcells == 2
            mutID = addmutations!(cellmodule.cells[parentcellid], cellmodule.cells[end], μ, 
                mutID, rng, mutationdist, Δt)
        else
            mutID = addmutations!(cellmodule.cells[parentcellid], μ, 
                mutID, rng, mutationdist, Δt)
        end
    end
    updatetime!(cellmodule, t)
    return cellmodule, subclones, mutID
end

function cellmutation!(cellmodule, subclones, selectioncoefficient, mutatingcell, t)
    
    #add new clone
    subcloneid = length(subclones) + 1
    parentid = getclonetype(mutatingcell)
    wildtype_rates = getwildtyperates(subclones)
    birthrate, deathrate, moranrate, asymmetricrate = get_newsubclone_rates(wildtype_rates, selectioncoefficient)
    newsubclone = Subclone(subcloneid, parentid, t, 1, birthrate, deathrate, moranrate, asymmetricrate)
    push!(subclones, newsubclone)

    #change clone type of new cell and update clone sizes
    setclonetype(mutatingcell, subcloneid)
    if parentid != 0
        subclones[parentid].size -= 1
    end
    return cellmodule, subclones
end

function celldeath!(cellmodule::CellModule, subclones, deadcell::Int64, args...)
    #frequency of cell type decreases
    clonetype = cellmodule.cells[deadcell].clonetype 
    subclones[clonetype].size -= 1
    #remove deleted cell
    deleteat!(cellmodule.cells, deadcell)

    return cellmodule
end

function celldeath!(cellmodule::CellModule, subclones, deadcells::Vector{Int64}, args...)
    for deadcell in deadcells
        clonetype = cellmodule.cells[deadcell].clonetype 
        subclones[clonetype].size -= 1
    end
    deleteat!(cellmodule.cells, sort(deadcells))

    return cellmodule
end

getnextID(population::SinglelevelPopulation) = getnextID(population.singlemodule.cells)

function getnextID(cells::CellVector)
    if all(no_mutations.(cells))
        return 1
    else
        allmutations = reduce(vcat, [cell.mutations for cell in cells])
        return maximum(allmutations)
    end
end

function getnextID(cells::AbstractTreeCellVector)
    nextID = 1
    for cellnode in cells
        if isnothing(cellnode) continue end
        if id(cellnode) + 1 > nextID
            nextID = id(cellnode) + 1
        end
    end
    return nextID
end

function copycell(cellold::Cell)
    return Cell(
        copy(cellold.mutations), 
        cellold.clonetype, 
        cellold.birthtime, 
        cellold.id, 
        cellold.parentid
    )
  end

function discretetime(rng, λ=1)
    return 1/λ
end

function exptime(rng::AbstractRNG)
    rand(rng, Exponential(1))
end

function exptime(rng::AbstractRNG, λ)
    rand(rng, Exponential(1/λ))
end

function no_mutations(cell)
    return length(cell.mutations) == 0
end

age(abstractmodule::AbstractModule) = abstractmodule.t
age(population::Population) = maximum(map(age, population))
age(population::SinglelevelPopulation) = age(population.singlemodule)
age(simulation::Simulation) = age(simulation.output)
age(multisim::MultiSimulation) = maximum(age(output) for output in multisim.output)

create_modulestructure(WellMixed, N) = WellMixed()
