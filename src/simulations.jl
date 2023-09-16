function simulate!(cellmodule::CellModule, input::BranchingMoranInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG; timefunc=exptime, t0=nothing, tmax=nothing)

    branchingprocess!(
        cellmodule, 
        input.birthrate, 
        input.deathrate, 
        input.Nmax, 
        input.μ, 
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        rng;
        timefunc,
        t0,
        numclones=input.numclones,
        mutant_selection=input.mutant_selection, 
        mutant_time=input.mutant_time, 
    )
    
    moranprocess!(
        cellmodule, 
        input.moranrate, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        input.μ, 
        input.mutationdist, 
        rng;
        timefunc,
        t0,
        numclones=input.numclones, 
        mutant_selection=input.mutant_selection, 
        mutant_time=input.mutant_time,
        moranincludeself=input.moranincludeself
    )

    return cellmodule
end

function simulate!(cellmodule::CellModule, input::BranchingInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG; timefunc=exptime, t0=nothing, tmax=nothing)

    branchingprocess!(
        cellmodule, 
        input.birthrate, 
        input.deathrate, 
        input.Nmax, 
        input.μ, 
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        rng; 
        timefunc,
        t0,
        numclones=input.numclones,
        mutant_selection=input.mutant_selection, 
        mutant_time=input.mutant_time, 
    )
    
    return cellmodule
end

function simulate!(cellmodule::CellModule, input::MoranInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG; timefunc=exptime, t0=nothing, tmax=nothing)


    moranprocess!(
        cellmodule, 
        input.moranrate, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        input.μ, 
        input.mutationdist, 
        rng; 
        timefunc,
        t0,
        numclones=input.numclones, 
        mutant_selection=input.mutant_selection, 
        mutant_time=input.mutant_time,
        moranincludeself=input.moranincludeself
    )

    return cellmodule
end

"""
    runsimulation_clonalmuts(::Type{Cell}, input::MoranInput[, rng::AbstractRNG]; tstep)

    Simulate a population of cells according to a Moran process for fixed time. Return 
    a vector of the number of clonal mutations acquired at given time intervals.
"""

function runsimulation_clonalmuts(::Type{Cell}, input::MoranInput, tstep, rng::AbstractRNG=Random.GLOBAL_RNG)

    cellmodule = initialize(Cell, input, rng)

    clonalmuts = Int64[]
    for t in 0:tstep:input.tmax
        cellmodule = moranprocess!(
            cellmodule, 
            input.moranrate, 
            input.tmax, 
            input.μ, 
            input.mutationdist, 
            rng; 
            numclones=input.numclones, 
            mutant_selection=input.mutant_selection, 
            mutant_time=input.mutant_time,
            moranincludeself=input.moranincludeself
        )
        push!(clonalmuts, clonal_mutations(cellmodule))
    end

    return Simulation(input,cellmodule), clonalmuts
end

"""
    branchingprocess!(cellmodule::CellModule, birthrate, deathrate, Nmax, μ, mutationdist, tmax,
    rng::AbstractRNG; numclones=0, mutant_selection=Float64[], mutant_time=Float64[])

Run branching process simulation, starting in state defined by cellmodule.

Simulate a stochastic branching process, starting with a single cell, with with birth rate 
`b`, death rate `d` until population reaches size `Nmax`.

Cells accumulate neutral mutations at division with rate `μ`.

If `numclones` = 0, all cells have the same fitness and there is only one (sub)clone. 
Otherwise, `numclones` is the number of fit subclones. The `i`th subclone arises by a single 
cell mutating at time `mutant_time[i]` and has selection coefficient `mutant_selection[i]`.

"""
function branchingprocess!(cellmodule::CellModule, birthrate, deathrate, Nmax, μ, mutationdist, tmax,
    rng::AbstractRNG; numclones=0, mutant_selection=Float64[], mutant_time=Float64[],timefunc=exptime, t0=nothing)
    
    t = !isnothing(t0) ? t0 : cellmodule.t
    N = length(cellmodule)
    
    mutID = N == 1 ? 1 : getnextID(cellmodule.cells)

    nclonescurrent = length(cellmodule.subclones) + 1  
    changemutrate = BitArray(undef, numclones + 1)
    changemutrate .= 1

    birthrates, deathrates = set_branching_birthdeath_rates( birthrate, deathrate, mutant_selection)

    #Rmax starts with birthrate + deathrate and changes once a fitter mutant is introduced, this ensures
    #that birthrate and deathrate have correct units
    Rmax = (maximum(birthrates[1:nclonescurrent])
                                + maximum(deathrates[1:nclonescurrent]))

    while N < Nmax

        #calc next event time and break if it exceeds tmax 
        Δt =  1/(Rmax * N) .* timefunc(rng)
        t = t + Δt
        if t > tmax
            break
        end

        randcell = rand(rng,1:N) #pick a random cell
        r = rand(rng,Uniform(0,Rmax))
        #get birth and death rates for randcell
        br = birthrates[cellmodule.cells[randcell].clonetype]
        dr = deathrates[cellmodule.cells[randcell].clonetype]

        if r < br 
            #cell divides
            cellmodule, mutID = celldivision!(cellmodule, randcell, t, mutID, μ, mutationdist, rng)
            N += 1
            #check if t>=mutant_time for next fit subclone
            if nclonescurrent < numclones + 1 && t >= mutant_time[nclonescurrent]
                #if current number clones != final number clones, one of the new cells is
                #mutated to become fitter and form a new clone
                    cellmodule, nclonescurrent = cellmutation!(cellmodule, N, N, t, 
                        nclonescurrent)
                
                    #change Rmax now there is a new fitter mutant
                    Rmax = (maximum(birthrates[1:nclonescurrent])
                                + maximum(deathrates[1:nclonescurrent]))
            end
            updatetime!(cellmodule, t)


        elseif r < br + dr
            #cell dies
            cellmodule = celldeath!(cellmodule, randcell)
            N -= 1
            updatetime!(cellmodule, t)
            #return empty cellmodule if all cells have died
            if N == 0
                return cellmodule
            end
        end
    end
    return cellmodule
end

"""
    getclonesize(cellmodule::CellModule)

Return number of cells in each subclone (including wild-type).
"""
function getclonesize(cellmodule::CellModule)
    return getclonesize(length(cellmodule), cellmodule.subclones)
end

"""
    getclonesize(N, subclones)
"""
function getclonesize(N, subclones)
   sizevec = [clone.size for clone in subclones]
   prepend!(sizevec, N - sum(sizevec)) 
end

"""
    set_branching_birthdeath_rates( birthrate, deathrate, mutant_selection)

Return Vectors of birthrates and deathrates for each subclone (including wild-type).
"""
function set_branching_birthdeath_rates(birthrate, deathrate, mutant_selection)
    birthrates = [birthrate]
    deathrates = [deathrate]
    #add birth and death rates for each subclone. 
    for i in 1:length(mutant_selection)
        push!(deathrates, deathrate)
        push!(birthrates, (1 + mutant_selection[i]) .* birthrate)
    end
    return birthrates,deathrates
end

function set_all_rates(birthrate, deathrate, mutant_selection)
    birthrates = [birthrate]
    deathrates = [deathrate]
    #add birth and death rates for each subclone. 
    for i in 1:length(mutant_selection)
        push!(deathrates, deathrate)
        push!(birthrates, (1 + mutant_selection[i]) .* birthrate)
    end
    return birthrates,deathrates
end

"""
    moranprocess!(cellmodule::CellModule, moranrate, tmax, μ, mutationdist, 
    rng::AbstractRNG; numclones=0, mutant_selection=Float64[], mutant_time=Float64[])

Run Moran process simulation, starting in state defined by cellmodule.

See also [`moranprocess`](@ref)

"""
function moranprocess!(cellmodule::CellModule, moranrate, tmax, μ, mutationdist, 
    rng::AbstractRNG; numclones=0, mutant_selection=Float64[], mutant_time=Float64[],
    timefunc=exptime, t0=nothing, moranincludeself=true)

    t = !isnothing(t0) ? t0 : cellmodule.t
    N = length(cellmodule)
    mutID = getnextID(cellmodule.cells)

    nclonescurrent = length(cellmodule.subclones) + 1  

    while true

        #calc next event time and break if it exceeds tmax 
        Δt =  1/(moranrate*N) .* timefunc(rng)
        t = t + Δt
        if t > tmax
            break
        end

        #pick a cell to divide proportional to clonetype
        if nclonescurrent == 1
            dividecellidx = rand(rng, 1:N)
        else
            p = [cell.clonetype==1 ? 1 : 1 + mutant_selection[cell.clonetype - 1] 
                    for cell in cellmodule.cells] 
            p /= sum(p)
            dividecellidx = sample(rng, 1:N, ProbabilityWeights(p)) 
        end

        #pick a random cell to die
        deadcellidx = if moranincludeself 
            deadcellidx = rand(rng, 1:N)
            #if dead cell and divide cell are the same kill one of the offspring
            deadcellidx = deadcellidx == dividecellidx ? N : deadcellidx
        else
            #exclude dividecellidx
            deadcellidx = rand(rng, deleteat!(collect(1:N), dividecellidx))
        end
        #cell divides
        cellmodule, mutID = celldivision!(cellmodule, dividecellidx, t, mutID, μ, mutationdist, rng)
        
        #check if t>=mutant_time for next fit subclone and more subclones are expected
        if nclonescurrent < numclones + 1 && t >= mutant_time[nclonescurrent]
            cellmodule, nclonescurrent = cellmutation!(cellmodule, N+1, N, t, nclonescurrent)
        end

        #cell dies
        cellmodule = celldeath!(cellmodule, deadcellidx)

        updatetime!(cellmodule, t)

    end
    return cellmodule
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
        return rand(rng,Poisson(μ))
    elseif mutationdist == :geometric
        return rand(rng, Geometric(1/(1+μ)))
    elseif mutationdist == :poissontimedep
        return rand(rng,Poisson(μ*Δt))
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

function celldivision!(cellmodule::CellModule, parentcell, t, mutID, μ, mutationdist, rng; nchildcells=2)
    
    Δt = t - cellmodule.cells[parentcell].birthtime
    cellmodule.cells[parentcell].birthtime = t
    if nchildcells == 2
        push!(cellmodule.cells, copycell(cellmodule.cells[parentcell])) #add new copy of parent cell to cells
        cellmodule.cells[end].id = cellmodule.cells[end-1].id + 1
        cellmodule.cells[end].parentid = cellmodule.cells[parentcell].id
    end
    #add new mutations to both new cells
    if μ > 0.0 
        if nchildcells == 2
            mutID = addmutations!(cellmodule.cells[parentcell], cellmodule.cells[end], μ, 
                mutID, rng, mutationdist, Δt)
        else
            mutID = addmutations!(cellmodule.cells[parentcell], μ, 
                mutID, rng, mutationdist, Δt)
        end
    end
    clonetype = cellmodule.cells[parentcell].clonetype
    if clonetype > 1 && nchildcells == 2
        cellmodule.subclones[clonetype - 1].size += 1
    end
    updatetime!(cellmodule, t)
    return cellmodule, mutID
end

function cellmutation!(cellmodule::CellModule, mutatingcell, N, t, nclonescurrent)
    
    #add new clone
    parenttype = cellmodule.cells[mutatingcell].clonetype
    mutations = deepcopy(cellmodule.cells[mutatingcell].mutations)
    Ndivisions = length(cellmodule.cells[mutatingcell].mutations)
    avdivisions = mean(map(x -> length(x.mutations), cellmodule.cells))
    push!(cellmodule.subclones, clone)

    #change clone type of new cell and update clone sizes
    nclonescurrent += 1
    cellmodule.cells[mutatingcell].clonetype = nclonescurrent

    if parenttype > 1
        cellmodule.subclones[parenttype - 1].size -= 1
    end


    return (;cellmodule, nclonescurrent)
end

function celldeath!(cellmodule::CellModule, deadcell::Int64, args...)
    #frequency of cell type decreases
    clonetype = cellmodule.cells[deadcell].clonetype 
    if clonetype > 1
        cellmodule.subclones[clonetype - 1].size -= 1
    end
    #remove deleted cell
    deleteat!(cellmodule.cells, deadcell)

    return cellmodule
end

function celldeath!(cellmodule::CellModule, deadcells::Vector{Int64}, args...)
    for deadcell in deadcells
        clonetype = cellmodule.cells[deadcell].clonetype 
        if clonetype > 1
            cellmodule.subclones[clonetype - 1].size -= 1
        end
    end
    deleteat!(cellmodule.cells, sort(deadcells))

    return cellmodule
end

cellremoval!(cellmodule::CellModule, deadcells::Vector{Int64}) = 
    celldeath!(cellmodule, deadcells)


function getnextID(cells::CellVector)
    if all(no_mutations.(cells))
        return 1
    else
        allmutations = reduce(vcat, [cell.mutations for cell in cells])
        return maximum(allmutations)
    end
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
age(population::Vector{T}) where T<:AbstractModule = maximum(map(age, population))
age(simulation::Simulation) = age(simulation.output)
age(multisim::MultiSimulation) = maximum(age(output) for output in multisim.output)

create_modulestructure(WellMixed, N) = WellMixed()
