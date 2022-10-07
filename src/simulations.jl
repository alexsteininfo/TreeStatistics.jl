"""
    runsimulation(input::SinglelevelInput, rng::AbstractRNG=Random.GLOBAL_RNG)
"""
function runsimulation(input::SinglelevelInput, rng::AbstractRNG=Random.GLOBAL_RNG)
    return runsimulation(Cell, input, rng)
end

function runsimulation(::Type{Cell}, input::SinglelevelInput, rng::AbstractRNG=Random.GLOBAL_RNG)
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later. 
    #UNLESS input.mutationdist=:poissontimedep or :fixedtimedep

    μ, clonalmutations, mutationdist = input.μ, input.clonalmutations, input.mutationdist
    if !(input.mutationdist ∈ (:poissontimedep, :fixedtimedep))
        input = newinput(input, μ=1, clonalmutations=0, mutationdist=:fixed)
    end

    cellmodule = initialize(Cell, input, rng)

    simulate!(cellmodule, input, rng)
    
    #If we set μ=1 etc we need to add proper mutations now (also remove undetectable subclones).
    if !(input.mutationdist ∈ (:poissontimedep, :fixedtimedep) )  
        processresults!(cellmodule, μ, clonalmutations, rng; mutationdist)
        input = newinput(input; μ, clonalmutations, mutationdist)
    #Otherwise if mutation accumulation is time dependent, add the final mutations.
    else
        final_timedep_mutations!(cellmodule::CellModule, input.μ, input.mutationdist, rng)
    end
    return Simulation(input, cellmodule)
end

function simulate!(cellmodule::CellModule, input::BranchingMoranInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG; timefunc=exptime, t0=nothing, tmax=nothing)

    branchingprocess!(
        cellmodule, 
        input.b, 
        input.d, 
        input.Nmax, 
        input.μ, 
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        rng;
        timefunc,
        t0,
        numclones=input.numclones,
        selection=input.selection, 
        tevent=input.tevent, 
        maxclonesize=Inf
    )
    
    moranprocess!(
        cellmodule, 
        input.bdrate, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        input.μ, 
        input.mutationdist, 
        rng;
        timefunc,
        t0,
        numclones=input.numclones, 
        selection=input.selection, 
        tevent=input.tevent
    )

    return cellmodule
end

function simulate!(cellmodule::CellModule, input::BranchingInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG; timefunc=exptime, t0=nothing, tmax=nothing)

    branchingprocess!(
        cellmodule, 
        input.b, 
        input.d, 
        input.Nmax, 
        input.μ, 
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        rng; 
        timefunc,
        t0,
        numclones=input.numclones,
        selection=input.selection, 
        tevent=input.tevent, 
        maxclonesize=input.maxclonesize
    )
    
    return cellmodule
end

function simulate!(cellmodule::CellModule, input::MoranInput, 
    rng::AbstractRNG=Random.GLOBAL_RNG; timefunc=exptime, t0=nothing, tmax=nothing)


    moranprocess!(
        cellmodule, 
        input.bdrate, 
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        input.μ, 
        input.mutationdist, 
        rng; 
        timefunc,
        t0,
        numclones=input.numclones, 
        selection=input.selection, 
        tevent=input.tevent
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
            input.bdrate, 
            input.tmax, 
            input.μ, 
            input.mutationdist, 
            rng; 
            numclones=input.numclones, 
            selection=input.selection, 
            tevent=input.tevent)
        push!(clonalmuts, clonal_mutations(cellmodule))
    end
    return cellmodule, clonalmuts


    return Simulation(input,cellmodule), clonalmuts
end

"""
    branchingprocess!(cellmodule::CellModule, b, d, Nmax, μ, mutationdist, tmax,
    rng::AbstractRNG; numclones=0, selection=Float64[], tevent=Float64[], maxclonesize=200)

Run branching process simulation, starting in state defined by cellmodule.

Simulate a stochastic branching process, starting with a single cell, with with birth rate 
`b`, death rate `d` until population reaches size `Nmax`.

Cells accumulate neutral mutations at division with rate `μ`, until all subclones exceed
`maxclonesize`.

If `numclones` = 0, all cells have the same fitness and there is only one (sub)clone. 
Otherwise, `numclones` is the number of fit subclones. The `i`th subclone arises by a single 
cell mutating at time `tevent[i]` and has selection coefficient `selection[i]`.

"""
function branchingprocess!(cellmodule::CellModule, b, d, Nmax, μ, mutationdist, tmax,
    rng::AbstractRNG; numclones=0, selection=Float64[], tevent=Float64[], maxclonesize=200,
    timefunc=exptime, t0=nothing)
    
    t = !isnothing(t0) ? t0 : cellmodule.tvec[end]
    N = cellmodule.Nvec[end]
    
    mutID = N == 1 ? 1 : getnextID(cellmodule.cells)

    nclonescurrent = length(cellmodule.subclones) + 1  
    executed = false
    changemutrate = BitArray(undef, numclones + 1)
    changemutrate .= 1

    birthrates, deathrates = set_branching_birthdeath_rates(b, d, selection)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures
    #that b and d have correct units
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
            #check if t>=tevent for next fit subclone
            if nclonescurrent < numclones + 1 && t >= tevent[nclonescurrent]
                #if current number clones != final number clones, one of the new cells is
                #mutated to become fitter and form a new clone
                    cellmodule, nclonescurrent = cellmutation!(cellmodule, N, N, t, 
                        nclonescurrent)
                
                    #change Rmax now there is a new fitter mutant
                    Rmax = (maximum(birthrates[1:nclonescurrent])
                                + maximum(deathrates[1:nclonescurrent]))
            end
            cellmodule = update_time_popsize!(cellmodule, t, N)


        elseif r < br + dr
            #cell dies
            cellmodule = celldeath!(cellmodule, randcell)
            N -= 1
            cellmodule = update_time_popsize!(cellmodule, t, N)
            #return empty cellmodule if all cells have died
            if N == 0
                return cellmodule
            end
        end

        #if population of all clones is sufficiently large no new mutations
        #are acquired, can use this approximation as only mutations above 1%
        #frequency can be reliably detected
        if ((maxclonesize !== nothing && 
            executed == false && 
            (getclonesize(cellmodule) .> maxclonesize) == changemutrate))

            μ = 0
            mutationdist=:fixed
            executed = true
        end
    end
    return cellmodule
end

"""
    getclonesize(cellmodule::CellModule)

Return number of cells in each subclone (including wild-type).
"""
function getclonesize(cellmodule::CellModule)
    return getclonesize(cellmodule.Nvec[end], cellmodule.subclones)
end

"""
    getclonesize(N, subclones)
"""
function getclonesize(N, subclones)
   sizevec = [clone.size for clone in subclones]
   prepend!(sizevec, N - sum(sizevec)) 
end

"""
    update_time_popsize(cellmodule::CellModule, t, N)

Update cellmodule with new time `t` and pop size `N`.
"""
function update_time_popsize!(cellmodule::CellModule, t, N)
    push!(cellmodule.tvec,t)
    push!(cellmodule.Nvec, N)
    return cellmodule
end

"""
    set_branching_birthdeath_rates(b, d, selection)

Return Vectors of birthrates and deathrates for each subclone (including wild-type).
"""
function set_branching_birthdeath_rates(b, d, selection)
    birthrates = [b]
    deathrates = [d]
    #add birth and death rates for each subclone. 
    for i in 1:length(selection)
        push!(deathrates, d)
        push!(birthrates, (1 + selection[i]) .* b)
    end
    return birthrates,deathrates
end


"""
    moranprocess!(cellmodule::CellModule, bdrate, tmax, μ, mutationdist, 
    rng::AbstractRNG; numclones=0, selection=Float64[], tevent=Float64[])

Run Moran process simulation, starting in state defined by cellmodule.

See also [`moranprocess`](@ref)

"""
function moranprocess!(cellmodule::CellModule, bdrate, tmax, μ, mutationdist, 
    rng::AbstractRNG; numclones=0, selection=Float64[], tevent=Float64[],
    timefunc=exptime, t0=nothing)

    t = !isnothing(t0) ? t0 : cellmodule.tvec[end]
    N = cellmodule.Nvec[end]
    mutID = getnextID(cellmodule.cells)

    nclonescurrent = length(cellmodule.subclones) + 1  

    while true

        #calc next event time and break if it exceeds tmax 
        Δt =  1/(bdrate*N) .* timefunc(rng)
        t = t + Δt
        if t > tmax
            break
        end

        #pick a cell to divide proportional to clonetype
        if nclonescurrent == 1
            randcelldivide = rand(rng, 1:N)
        else
            p = [cell.clonetype==1 ? 1 : 1 + selection[cell.clonetype - 1] 
                    for cell in cellmodule.cells] 
            p /= sum(p)
            randcelldivide = sample(rng, 1:N, ProbabilityWeights(p)) 
        end

        #pick a random cell to die
        randcelldie = rand(rng,1:N) 

        #cell divides
        cellmodule, mutID = celldivision!(cellmodule, randcelldivide, t, mutID, μ, mutationdist, rng)
        
        #check if t>=tevent for next fit subclone and more subclones are expected
        if nclonescurrent < numclones + 1 && t >= tevent[nclonescurrent]
            cellmodule, nclonescurrent = cellmutation!(cellmodule, N+1, N, t, nclonescurrent)
        end

        #cell dies
        cellmodule = celldeath!(cellmodule, randcelldie)

        cellmodule = update_time_popsize!(cellmodule, t, N)

    end
    return cellmodule
end

"""
    initialize(::Type{T}, input, rng::AbstractRNG=Random.GLOBAL_RNG)

Initialise population of cells based on `input` and return as a `CellModule` if 
`T` is a `Cell` or a `TreeModule` if `T` is an `AbstractTreeCell`.
"""
function initialize end

function initialize(::Type{Cell}, input::Union{BranchingInput, BranchingMoranInput}, 
    rng::AbstractRNG=Random.GLOBAL_RNG)
    
    return initializesim_branching(
        input.Nmax,
        clonalmutations=input.clonalmutations,
    )
end


function initialize(::Type{Cell}, input::MoranInput, rng::AbstractRNG=Random.GLOBAL_RNG)

    return initializesim_moran(
        input.N, 
        clonalmutations=input.clonalmutations,
    )
end


function initializesim_branching(Nmax=nothing; clonalmutations=0, id=1, parentid=0)

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #population starts with one cell
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)

    #Initialize array of cell type that stores mutations for each cell and their clone type
    #clone type of 1 is the host population with selection=0
    cells = Cell[]
    if Nmax !== nothing 
        sizehint!(cells, Nmax)
    end
    push!(cells, Cell([], 1, 0, id, parentid))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    mutID = addnewmutations!(cells[1], clonalmutations, mutID)

    subclones = CloneTracker[]

    cellmodule = CellModule(
        Nvec,
        tvec,
        cells,
        subclones,
        1,
        0
    )
    return cellmodule
end

function initializesim_moran(N; clonalmutations=0)

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #Initialize array of cell type that stores mutations for each cell and their clone type
    #clone type of 1 is the host population with selection=0
    cells = [Cell(Int64[], 1, 0, id, 0) for id in 1:N]

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    for cell in cells
        cell.mutations = collect(1:clonalmutations)
    end

    subclones = CloneTracker[]

    cellmodule = CellModule(
        [N],
        tvec,
        cells,
        subclones,
        1,
        0
    )
    return cellmodule
end

function initialize_from_cells(::Type{T}, cells, subclones::Vector{CloneTracker}, 
    id, parentid; inittime=0.0) where T <: AbstractModule

    #population starts from list of cells
    tvec = Float64[inittime]
    Nvec = Int64[length(cells)]

    cellmodule = T(
        Nvec,
        tvec,
        cells,
        subclones,
        id,
        parentid
    )
    return cellmodule
end

function initialize_from_cells(cells::Vector{Cell}, subclones::Vector{CloneTracker}, id, parentid; inittime=0.0)
    return initialize_from_cells(
        CellModule, 
        cells, 
        subclones,
        id, 
        parentid; 
        inittime
    )
end

function addmutations!(cell1, cell2, μ, mutID, rng, mutationdist=mutationdist, Δt=Δt)
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

# function newmutations!(cell, μ, mutID, rng::AbstractRNG; mutationdist=:poisson, Δt=nothing)
#     #function to add new mutations to cells based on μ
#     numbermutations = numbernewmutations(rng, mutationdist, μ, Δt=Δt)
#     return addnewmutations!(cell, numbermutations, mutID)
# end

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

function celldivision!(cellmodule::CellModule, parentcell, t, mutID, μ, mutationdist, rng)
    
    Δt = t - cellmodule.cells[parentcell].birthtime
    cellmodule.cells[parentcell].birthtime = t
    push!(cellmodule.cells, copycell(cellmodule.cells[parentcell])) #add new copy of parent cell to cells
    cellmodule.cells[end].id = cellmodule.cells[end-1].id + 1
    cellmodule.cells[end].parentid = cellmodule.cells[parentcell].id
    #add new mutations to both new cells
    if μ > 0.0 
        mutID = addmutations!(cellmodule.cells[parentcell], cellmodule.cells[end], μ, 
            mutID, rng, mutationdist, Δt)
    end
    clonetype = cellmodule.cells[parentcell].clonetype
    if clonetype > 1
        cellmodule.subclones[clonetype - 1].size += 1
    end

    return cellmodule, mutID
end

function cellmutation!(cellmodule::CellModule, mutatingcell, N, t, nclonescurrent)
    
    #add new clone
    parenttype = cellmodule.cells[mutatingcell].clonetype
    mutations = deepcopy(cellmodule.cells[mutatingcell].mutations)
    Ndivisions = length(cellmodule.cells[mutatingcell].mutations)
    avdivisions = mean(map(x -> length(x.mutations), cellmodule.cells))
    clone = CloneTracker(parenttype, cellmodule.id, t, mutations, N, Ndivisions, 
        avdivisions, 1)
    push!(cellmodule.subclones, clone)

    #change clone type of new cell and update clone sizes
    nclonescurrent += 1
    cellmodule.cells[mutatingcell].clonetype = nclonescurrent

    if parenttype > 1
        cellmodule.subclones[parenttype - 1].size -= 1
    end


    return cellmodule, nclonescurrent
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


function getnextID(cells::Vector{Cell})
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

age(abstractmodule::AbstractModule) = abstractmodule.tvec[end]
age(population::Vector{T}) where T<:AbstractModule = maximum(map(age, population))
age(simulation::Simulation) = age(simulation.output)
age(multisim::MultiSimulation) = maximum(age(output) for output in multisim.output)