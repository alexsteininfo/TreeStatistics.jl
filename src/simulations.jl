"""
    simulate!(population, input, selection, counters, rng::AbstractRNG=Random.GLOBAL_RNG;
        timefunc=exptime, t0=age(population), tmax=input.tmax)

Run a simulation defined by `input` and `selection`, with `population` giving the
initial state. The time between update events is determined by `timefunc`. Simulation time
begins at `t0` and runs to a maximum time given by `minimum((tmax, input.tmax))`.

"""
function simulate! end

function simulate!(
    population::SinglelevelPopulation,
    input::BranchingMoranInput,
    selection::AbstractSelection,
    nextID::Integer,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    t0=nothing,
    tmax=nothing
)
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
    #add final time-dependent mutations
    final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=input.tmax)
    return population, nextID
end

function simulate!(
    population::SinglelevelPopulation,
    input::BranchingInput,
    selection::AbstractSelection,
    nextID::Integer,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    t0=nothing,
    tmax=nothing
)
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
    #add final time-dependent mutations
    final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=input.tmax)
    return population, nextID
end

function simulate!(
    population::SinglelevelPopulation,
    input::MoranInput,
    selection::AbstractSelection,
    nextID::Integer,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    t0=nothing,
    tmax=nothing
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
    #add final time-dependent mutations
    final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=input.tmax)
    return population, nextID
end

"""
    branchingprocess!(
        population::SinglelevelPopulation,
        selection::AbstractSelection,
        Nmax,
        μ,
        mutationdist,
        tmax,
        nextID,
        rng::AbstractRNG;
        timefunc=exptime,
        t0=nothing
    )

Run branching process simulation, starting in state defined by cellmodule.

Simulate a stochastic branching process, starting with a single cell, with with birth rate
`b`, death rate `d` until population reaches size `Nmax`.

Cells accumulate neutral mutations at division with rate `μ`.
"""
function branchingprocess!(
    population::SinglelevelPopulation,
    selection::AbstractSelection,
    Nmax,
    μ,
    mutationdist,
    tmax,
    nextID,
    rng::AbstractRNG;
    timefunc=exptime,
    t0=nothing
)
    t = !isnothing(t0) ? t0 : age(population.singlemodule)
    N = length(population.singlemodule)

    nsubclones = getmaxsubclones(selection)
    nsubclonescurrent = length(population.subclones)
    birthrates = getbirthrates(population.subclones)
    deathrates = getdeathrates(population.subclones)

    #Rmax starts with birthrate + deathrate and changes once a fitter mutant is introduced
    Rmax = maximum(birthrates) + maximum(deathrates)

    if(birthrate>deathrate)
        while N < Nmax && N > 0

            #calc next event time and break if it exceeds tmax
            Δt =  1 / (Rmax * N) .* timefunc(rng)
            t + Δt <= tmax || break # end simulation if time exceeds maximum
            t += Δt

            population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones =
                branchingupdate!(
                    population,
                    selection,
                    birthrates,
                    deathrates,
                    Rmax,
                    N,
                    nextID,
                    nsubclonescurrent,
                    nsubclones,
                    t,
                    μ,
                    mutationdist,
                    rng
                )
        end
    elseif(birthrate<deathrate)
        while N > Nmax && N > 0

            #calc next event time and break if it exceeds tmax
            Δt =  1 / (Rmax * N) .* timefunc(rng)
            t + Δt <= tmax || break # end simulation if time exceeds maximum
            t += Δt

            population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones =
                branchingupdate!(
                    population,
                    selection,
                    birthrates,
                    deathrates,
                    Rmax,
                    N,
                    nextID,
                    nsubclonescurrent,
                    nsubclones,
                    t,
                    μ,
                    mutationdist,
                    rng
                )
        end
    else
        while N < Nmax && N > 0

            #calc next event time and break if it exceeds tmax
            Δt =  1 / (Rmax * N) .* timefunc(rng)
            t + Δt <= tmax || break # end simulation if time exceeds maximum
            t += Δt

            population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones =
                branchingupdate!(
                    population,
                    selection,
                    birthrates,
                    deathrates,
                    Rmax,
                    N,
                    nextID,
                    nsubclonescurrent,
                    nsubclones,
                    t,
                    μ,
                    mutationdist,
                    rng
                )
        end
    end
    return population, nextID
end

function branchingupdate!(
    population::SinglelevelPopulation,
    selection,
    birthrates,
    deathrates,
    Rmax,
    N,
    nextID,
    nsubclonescurrent,
    nsubclones,
    t,
    μ,
    mutationdist,
    rng
)
    randcellid = rand(rng, 1:N) #pick a random cell
    randcell = population.singlemodule.cells[randcellid]
    r = rand(rng, Uniform(0, Rmax))
    #get birth and death rates for randcell
    cellsubclone = getclonetype(randcell)
    br = birthrates[cellsubclone]
    dr = deathrates[cellsubclone]

    if r < br
        #cell divides
        cellmodule, subclones, nextID = celldivision!(
            population.singlemodule,
            population.subclones,
            randcellid,
            t,
            nextID,
            μ,
            mutationdist,
            rng
        )
        N += 1
        #check if ready for a new mutant subclone
        if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
            newmutant_selectioncoeff =
                getselectioncoefficient(selection, nsubclonescurrent, rng)
            cellmutation!(
                cellmodule,
                population.subclones,
                newmutant_selectioncoeff,
                cellmodule[randcellid],
                t
            )
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
    moranprocess!(
        population::SinglelevelPopulation,
        selection,
        μ,
        mutationdist,
        tmax,
        nextID,
        rng::AbstractRNG;
        timefunc=exptime,
        t0=nothing,
        moranincludeself=true
    )

Run Moran process simulation, starting in state defined by cellmodule.
"""
function moranprocess!(
    population::SinglelevelPopulation,
    selection,
    μ,
    mutationdist,
    tmax,
    nextID,
    rng::AbstractRNG;
    timefunc=exptime,
    t0=nothing,
    moranincludeself=true
)
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
        population, moranrates, Rmax, N, nextID, nsubclonescurrent = moranupdate!(
            population,
            selection,
            moranrates,
            Rmax,
            N,
            nextID,
            nsubclonescurrent,
            nsubclones,
            t,
            μ,
            mutationdist,
            moranincludeself,
            rng
        )
    end
    return population, nextID
end

function moranupdate!(
    population::SinglelevelPopulation,
    selection,
    moranrates,
    Rmax,
    N,
    nextID,
    nsubclonescurrent,
    nsubclones,
    t,
    μ,
    mutationdist,
    moranincludeself,
    rng
)
    dividecellid = rand(rng, 1:N) #pick a random cell to divide
    dividecell = population.singlemodule.cells[dividecellid]
    r = rand(rng, Uniform(0, Rmax))
    mr = moranrates[getclonetype(dividecell)]
    if r < mr
        #pick a random cell to die
        deadcell = choose_moran_deadcell(N, dividecellid, moranincludeself, rng)
        #cell divides
        cellmodule, subclones, nextID = celldivision!(
            population.singlemodule,
            population.subclones,
            dividecellid,
            t,
            nextID,
            μ,
            mutationdist,
            rng
        )
        #check if t>=mutant_time for next fit subclone and more subclones are expected
            if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
                newmutant_selectioncoeff =
                    getselectioncoefficient(selection, nsubclonescurrent, rng)
                cellmutation!(
                    cellmodule,
                    population.subclones,
                    newmutant_selectioncoeff,
                    population.singlemodule.cells[dividecellid],
                    t
                )
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


updatetime!(abstractmodule, t) = abstractmodule.t = t

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
age(population::MultilevelPopulation) = maximum(map(age, population))
age(population::SinglelevelPopulation) = age(population.singlemodule)
age(simulation::Simulation) = age(simulation.output)
age(multisim::MultiSimulation) = maximum(age(output) for output in multisim.output)
