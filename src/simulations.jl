"""
    run1simulation(input::BranchingMoranInput[, rng::AbstractRNG])

    Simulate a population of cells that grows by a branching process to fixed size then 
    switches to a Moran process.
    
    Take simulation parameters from `input` and return a Simulation.
"""
function run1simulation(input::BranchingMoranInput, rng::AbstractRNG = Random.GLOBAL_RNG)

    #Run branching simulation starting with a single cell.
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later.
    moduletracker = 
        branchingprocess(input.b, input.d, input.Nmax, 1, rng, numclones = input.numclones, 
            fixedmu = true, clonalmutations = 0, selection = input.selection,
            tevent = input.tevent, maxclonesize = Inf)
    
    moduletracker = 
        moranprocess!(moduletracker, input.bdrate, input.tmax, 1, rng::AbstractRNG; 
        numclones = input.numclones, fixedmu = true, selection = input.selection, 
        tevent = input.tevent)

    moduletracker = 
        processresults!(moduletracker, input.μ, input.clonalmutations, rng)
    
    return Simulation(input, moduletracker)
end

"""
    run1simulation(input::BranchingInput[, rng::AbstractRNG])

    Simulate a population of cells that grows by a branching process to a fixed size.
"""
function run1simulation(input::BranchingInput, rng::AbstractRNG = Random.GLOBAL_RNG)

    #Run branching simulation starting with a single cell.
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later.
    moduletracker = 
        branchingprocess(input.b, input.d, input.Nmax, 1, rng, numclones = input.numclones, 
            fixedmu = true, clonalmutations = 0, selection = input.selection,
            tevent = input.tevent, maxclonesize = input.maxclonesize)
    
    #Add mutations and process simulation output to get SimResults.
    #Remove undetectable subclones from moduletracker
    moduletracker = 
        processresults!(moduletracker, input.μ, input.clonalmutations, rng)

    return Simulation(input, moduletracker)
end

"""
    run1simulation(input::MoranInput[, rng::AbstractRNG])

    Simulate a population of cells according to a Moran process for fixed time.
"""

function run1simulation(input::MoranInput, rng::AbstractRNG = Random.GLOBAL_RNG)

    #Run Moran simulation starting from a population of N identical cells.
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later.
    moduletracker = 
        moranprocess(input.N, input.bdrate, input.tmax, 1, rng, 
                    numclones = input.numclones, fixedmu = true, 
                    clonalmutations = 0, selection = input.selection,
                    tevent = input.tevent)

    #Add mutations and process simulation output to get SimResults.
    #Remove undetectable subclones from moduletracker   
    moduletracker = 
        processresults!(moduletracker, input.μ, input.clonalmutations, rng)
    return Simulation(input,moduletracker)
end

"""
    branchingprocess(input::BranchingInput, rng::AbstractRNG; <keyword arguments>)

Simulate a stochastic branching process with parameters defined by input and return 
ModuleTracker. 

Simulation is by a rejection-kinetic Monte Carlo algorithm and starts with a single cell. 

"""
function branchingprocess(input::BranchingInput, rng::AbstractRNG,
                            fixedmu=input.fixedmu, μ=input.μ, 
                            clonalmutations=input.clonalmutations)

    return branchingprocess(input.b, input.d, input.Nmax, μ, rng, 
                            numclones = input.numclones, fixedmu = fixedmu, 
                            clonalmutations = clonalmutations, 
                            selection = input.selection, tevent = input.tevent, 
                            maxclonesize = input.maxclonesize)
end

"""
    branchingprocess(b, d, Nmax, μ, rng::AbstractRNG; <keyword arguments>)

Simulate a stochastic branching process, starting with a single cell, with with birth rate 
`b`, death rate `d` until population reaches size `Nmax`.

Cells accumulate neutral mutations at division with rate `μ`, until all subclones exceed
`maxclonesize`.

If `numclones` = 0, all cells have the same fitness and there is only one (sub)clone. 
Otherwise, `numclones` is the number of fit subclones. The `i`th subclone arises by a single 
cell mutating at time `tevent[i]` and has selection coefficient `selection[i]`.

"""
function branchingprocess(b, d, Nmax, μ, rng::AbstractRNG; numclones=0, fixedmu=false,
    clonalmutations=μ, selection=Float64[], tevent=Float64[], maxclonesize=200)

    #initialize arrays and parameters
    moduletracker = initializesim_branching(
        Nmax, 
        clonalmutations=clonalmutations
    )
    
    #run simulation
    moduletracker = branchingprocess!(moduletracker, b, d, Nmax, μ, rng, numclones=numclones,
        fixedmu=fixedmu, selection=selection, tevent=tevent, maxclonesize=maxclonesize)
    return moduletracker
end

"""
    branchingprocess!(moduletracker::ModuleTracker, b, d, Nmax, μ, rng::AbstractRNG; 
        <keyword arguments>)

Run branching process simulation, starting in state defined by moduletracker.

See also [`branchingprocess`](@ref)
"""
function branchingprocess!(moduletracker::ModuleTracker, b, d, Nmax, μ, rng::AbstractRNG; 
    numclones=0, fixedmu=false, selection=selection, tevent=Float64[], 
    maxclonesize=200, tmax=Inf)

    t, N = moduletracker.tvec[end], moduletracker.Nvec[end]
    mutID = N == 1 ? 1 : getmutID(moduletracker.cells)

    nclonescurrent = length(moduletracker.subclones) + 1  
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
        Δt =  1/(Rmax * N) .* exptime(rng)
        t = t + Δt
        if t > tmax
            break
        end

        randcell = rand(rng,1:N) #pick a random cell
        r = rand(rng,Uniform(0,Rmax))
        #get birth and death rates for randcell
        br = birthrates[moduletracker.cells[randcell].clonetype]
        dr = deathrates[moduletracker.cells[randcell].clonetype]

        if r < br 
            #cell divides
            moduletracker, mutID = celldivision!(moduletracker, randcell, mutID, μ, rng,
                fixedmu = fixedmu)
            N += 1
            #check if t>=tevent for next fit subclone
            if nclonescurrent < numclones + 1 && t >= tevent[nclonescurrent]
                #if current number clones != final number clones, one of the new cells is
                #mutated to become fitter and form a new clone
                    moduletracker, nclonescurrent = cellmutation!(moduletracker, N, N, t, 
                        nclonescurrent)
                
                    #change Rmax now there is a new fitter mutant
                    Rmax = (maximum(birthrates[1:nclonescurrent])
                                + maximum(deathrates[1:nclonescurrent]))
            end

        elseif r < br + dr
            #cell dies
            moduletracker = celldeath!(moduletracker, randcell)
            N -= 1
        end

        moduletracker = update_time_popsize!(moduletracker, t, N)

        #if population of all clones is sufficiently large no new mutations
        #are acquired, can use this approximation as only mutations above 1%
        #frequency can be reliably detected
        if ((maxclonesize !== nothing && 
            executed == false && 
            (getclonesize(moduletracker) .> maxclonesize) == changemutrate))

            μ = 0
            fixedmu = true
            executed = true
        end
    end
    return moduletracker
end

"""
    getclonesize(moduletracker::ModuleTracker)

Return number of cells in each subclone (including wild-type).
"""
function getclonesize(moduletracker::ModuleTracker)
    return getclonesize(moduletracker.Nvec[end], moduletracker.subclones)
end

"""
    getclonesize(N, subclones)
"""
function getclonesize(N, subclones)
   sizevec = [clone.size for clone in subclones]
   prepend!(sizevec, N - sum(sizevec)) 
end

"""
    update_time_popsize(moduletracker::ModuleTracker, t, N)

Update moduletracker with new time `t` and pop size `N`.
"""
function update_time_popsize!(moduletracker::ModuleTracker, t, N)
    push!(moduletracker.tvec,t)
    push!(moduletracker.Nvec, N)
    return moduletracker
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
    moranprocess(input::MoranInput, rng::AbstractRNG; <keyword arguments>)

Simulate a Moran process with parameters defined by input and return ModuleTracker. 

Simulation is by a Gillespie algorithm.

"""
function moranprocess(input::MoranInput, rng::AbstractRNG,
                    fixedmu=input.fixedmu, μ=input.μ, 
                    clonalmutations=input.clonalmutations)

    return moranprocess(input.N, input.bdrate, input.tmax, μ, rng, 
                        numclones = input.numclones, fixedmu = fixedmu, 
                        clonalmutations = clonalmutations, 
                        selection = input.selection, tevent = input.tevent)
end

"""
    moranprocess(N, bdrate, tmax, μ, rng::AbstractRNG; <keyword arguments>)

Simulate a Moran process starting with `N` cells until time `tmax`.

Update events comprise of a birth and a death, and occur with rate `bdrate`. Cells 
accumulate neutral mutations at division with rate `μ`, until all subclones exceed
`maxclonesize`.

If `numclones` = 0, all cells have the same fitness and there is only one (sub)clone. 
Otherwise, `numclones` is the number of fit subclones. The `i`th subclone arises by a single 
cell mutating at time `tevent[i]` and has selection coefficient `selection[i]`. 

"""
function moranprocess(N, bdrate, tmax, μ, rng::AbstractRNG; numclones = 0, fixedmu = false,
    clonalmutations = μ, selection = Float64[], tevent = Float64[])

    moduletracker = initializesim_moran(N, clonalmutations = clonalmutations)

    #run simulation
    moduletracker = moranprocess!(moduletracker, bdrate, tmax, μ, rng, numclones = numclones,
    fixedmu = fixedmu, selection = selection, tevent = tevent)
    return moduletracker
end
"""
    moranprocess!(moduletracker::ModuleTracker, bdrate, tmax, μ, rng::AbstractRNG; 
        <keyword-arguments>)

Run Moran process simulation, starting in state defined by moduletracker.

See also [`moranprocess`](@ref)

"""
function moranprocess!(moduletracker::ModuleTracker, bdrate, tmax, μ, rng::AbstractRNG; 
    numclones = 0, fixedmu = false, selection = Float64[], tevent = Float64[])

    t, N = moduletracker.tvec[end], moduletracker.Nvec[end]
    mutID = getmutID(moduletracker.cells)

    nclonescurrent = length(moduletracker.subclones) + 1  

    while true

        #calc next event time and break if it exceeds tmax 
        Δt =  1/(bdrate*N) .* exptime(rng)
        t = t + Δt
        if t > tmax
            break
        end

        #pick a cell to divide proportional to clonetype
        if nclonescurrent == 1
            randcelldivide = rand(rng, 1:N)
        else
            p = [cell.clonetype==1 ? 1 : 1 + selection[cell.clonetype - 1] 
                    for cell in moduletracker.cells] 
            p /= sum(p)
            randcelldivide = sample(rng, 1:N, ProbabilityWeights(p)) 
        end

        #pick a random cell to die
        randcelldie = rand(rng,1:N) 

        #cell divides
        moduletracker, mutID = celldivision!(moduletracker, randcelldivide, mutID, μ, rng,
            fixedmu = fixedmu)
        
        #check if t>=tevent for next fit subclone and more subclones are expected
        if nclonescurrent < numclones + 1 && t >= tevent[nclonescurrent]
            moduletracker, nclonescurrent = cellmutation!(moduletracker, N+1, N, t, nclonescurrent)
        end

        #cell dies
        moduletracker = celldeath!(moduletracker, randcelldie)

        moduletracker = update_time_popsize!(moduletracker, t, N)

    end
    return moduletracker
end

"""
    initializesim(input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG)

Initialise simulation and return a ModuleTracker.
"""
function initializesim(siminput::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG)
    
    return initializesim_branching(
        siminput.Nmax,
        clonalmutations=siminput.clonalmutations,
    )
end

"""
    initializesim(input::MoranInput, rng::AbstractRNG=Random.GLOBAL_RNG)
"""
function initializesim(siminput::MoranInput, rng::AbstractRNG=Random.GLOBAL_RNG)

    return initializesim_moran(
        siminput.N, 
        clonalmutations=siminput.clonalmutations,
    )
end

"""
    initializesim_branching(input::BranchingInput, rng::AbstractRNG=Random.GLOBAL_RNG)

Initialise simulation and return a ModuleTracker.
"""
function initializesim_branching(Nmax=nothing; clonalmutations=0)

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
    push!(cells, Cell([],1))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    cells[1],mutID = newmutations!(cells[1], clonalmutations, mutID)

    subclones = CloneTracker[]

    moduletracker = ModuleTracker(
        Nvec,
        tvec,
        cells,
        subclones,
        1,
        0
    )
    return moduletracker
end

function initializesim_moran(N; clonalmutations=0)

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #Initialize array of cell type that stores mutations for each cell and their clone type
    #clone type of 1 is the host population with selection=0
    cells = [Cell(Int64[],1) for _ in 1:N]

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    for cell in cells
        cell.mutations = collect(1:clonalmutations)
    end

    subclones = CloneTracker[]

    moduletracker = ModuleTracker(
        [N],
        tvec,
        cells,
        subclones,
        1,
        0
    )
    return moduletracker
end

function initializesim_from_cells(cells::Array{Cell,1}, subclones::Array{CloneTracker, 1}, 
    id, parentid;
    inittime=0.0)

    #initialize time to zero
    t = inittime
    tvec = Float64[]
    push!(tvec,t)

    #population starts from list of cells
    N = length(cells)
    Nvec = Int64[]
    push!(Nvec, N)

    for (i, subclone) in enumerate(subclones)
        subclone.size = sum(cell.clonetype .== 1)
    end

    moduletracker = ModuleTracker(
        Nvec,
        tvec,
        cells,
        subclones,
        id,
        parentid
    )
    return moduletracker
end


function newmutations!(cell, μ, mutID, rng::AbstractRNG; fixedmu = false)
    #function to add new mutations to cells based on μ
    numbermutations = fixedmu ? μ : rand(rng,Poisson(μ))
    return newmutations!(cell, numbermutations, mutID)
end

function newmutations!(cell, μ, mutID)
    #function to add new mutations to cells based on μ
    append!(cell.mutations, mutID:mutID + μ - 1)
    mutID = mutID + μ
    return cell, mutID
end

function celldivision!(moduletracker::ModuleTracker, parentcell, mutID, μ, 
    rng::AbstractRNG; fixedmu=true)
    
    push!(moduletracker.cells, copycell(moduletracker.cells[parentcell])) #add new copy of parent cell to cells
    #add new mutations to both new cells
    if μ > 0.0 
        moduletracker.cells[parentcell],mutID = newmutations!(moduletracker.cells[parentcell], μ, 
            mutID, rng, fixedmu = fixedmu)
        moduletracker.cells[end],mutID = newmutations!(moduletracker.cells[end], μ, mutID, rng,
            fixedmu=fixedmu)
    end
    clonetype = moduletracker.cells[parentcell].clonetype
    if clonetype > 1
        moduletracker.subclones[clonetype - 1].size += 1
    end

    return moduletracker, mutID
end

function cellmutation!(moduletracker::ModuleTracker, mutatingcell, N, t, nclonescurrent)
    
    #add new clone
    parenttype = moduletracker.cells[mutatingcell].clonetype
    mutations = deepcopy(moduletracker.cells[mutatingcell].mutations)
    Ndivisions = length(moduletracker.cells[mutatingcell].mutations)
    avdivisions = mean(map(x -> length(x.mutations), moduletracker.cells))
    clone = CloneTracker(parenttype, moduletracker.id, t, mutations, N, Ndivisions, 
        avdivisions, 1)
    push!(moduletracker.subclones, clone)

    #change clone type of new cell and update clone sizes
    nclonescurrent += 1
    moduletracker.cells[mutatingcell].clonetype = nclonescurrent

    if parenttype > 1
        moduletracker.subclones[parenttype - 1].size -= 1
    end


    return moduletracker, nclonescurrent
end

function celldeath!(moduletracker::ModuleTracker, deadcell::Int64)
    #frequency of cell type decreases
    clonetype = moduletracker.cells[deadcell].clonetype 
    if clonetype > 1
        moduletracker.subclones[clonetype - 1].size -= 1
    end
    #remove deleted cell
    deleteat!(moduletracker.cells, deadcell)

    return moduletracker
end

function celldeath!(moduletracker::ModuleTracker, deadcells::Array{Int64, 1})
    for deadcell in deadcells
        clonetype = moduletracker.cells[deadcell].clonetype 
        if clonetype > 1
            moduletracker.subclones[clonetype - 1].size -= 1
        end
    end
    deleteat!(moduletracker.cells, sort(deadcells))

    return moduletracker
end

function getmutID(cells::Vector{Cell})
    if all(no_mutations.(cells))
        return 1
    else
        allmutations = reduce(vcat, [cell.mutations for cell in cells])
        return maximum(allmutations)
    end
end

function copycell(cellold::Cell)
    return Cell(copy(cellold.mutations), cellold.clonetype)
  end

function exptime(rng::AbstractRNG)
    - log(rand(rng))
end

function no_mutations(cell)
    return length(cell.mutations) == 0
end
