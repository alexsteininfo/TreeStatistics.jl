"""
    run1simulation(IP::InputParameters{BranchingInput}, rng::AbstractRNG = MersenneTwister; 
        <keyword arguments>)

    Run a single branching process simulation using parameters defined by IP and return
    a Simulation object.
"""

function run1simulation(IP::InputParameters{BranchingInput}, rng::AbstractRNG = Random.GLOBAL_RNG;
    minclonefreq = 0.0, maxclonefreq = 1.0)

    #Run branching simulation starting with a single cell.
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later.
    simtracker = 
        branchingprocess(IP.siminput.b, IP.siminput.d, IP.siminput.Nmax, 1, rng, numclones = IP.siminput.numclones, 
            fixedmu = true, clonalmutations = 0, selection = IP.siminput.selection,
            tevent = IP.siminput.tevent, maxclonesize = IP.siminput.maxclonesize, 
            timefunction = IP.siminput.timefunction)
    
    #Add mutations and process simulation output to get SimResults.
    #Remove undetectable subclones from simtracker
    simtracker, simresults = 
        processresults!(simtracker, IP.siminput.Nmax, IP.siminput.numclones, IP.siminput.μ, 
                        IP.siminput.fixedmu, IP.siminput.clonalmutations, IP.ploidy, 
                        minclonefreq, maxclonefreq, rng)
    
    #Mimic experimental data by sampling from the true VAF
    sampleddata = 
        sampledhist(simresults.trueVAF, IP.siminput.Nmax, rng, 
                    detectionlimit = IP.detectionlimit, 
                    read_depth = IP.read_depth, cellularity = IP.cellularity)

    return Simulation(IP,simresults,sampleddata)
end

"""
    run1simulation(IP::InputParameters{MorangInput}, rng::AbstractRNG = MersenneTwister; 
        <keyword arguments>)

    Run a single moran process simulation using parameters defined by IP and return
    a Simulation object.
"""

function run1simulation(IP::InputParameters{MoranInput}, rng::AbstractRNG = Random.GLOBAL_RNG;
    minclonefreq = 0.0, maxclonefreq = 1.0)

    #Run Moran simulation starting from a population of N identical cells.
    #Initially set clonalmutations = 0 and μ = 1. These are expanded later.
    simtracker = 
        moranprocess(IP.siminput.N, IP.siminput.bdrate, IP.siminput.tmax, 1, rng, 
                    numclones = IP.siminput.numclones, fixedmu = true, 
                    clonalmutations = 0, selection = IP.siminput.selection,
                    tevent = IP.siminput.tevent, 
                    timefunction = IP.siminput.timefunction)

    #Add mutations and process simulation output to get SimResults.
    #Remove undetectable subclones from simtracker   
    simtracker, simresults = 
        processresults!(simtracker, IP.siminput.N, IP.siminput.numclones, IP.siminput.μ, 
                        IP.siminput.fixedmu, IP.siminput.clonalmutations, 
                        IP.ploidy, minclonefreq, maxclonefreq, rng)
    
    #Mimic experimental data by sampling from the true VAF
    sampleddata = 
        sampledhist(simresults.trueVAF, IP.siminput.N, rng, 
                    detectionlimit = IP.detectionlimit, 
                    read_depth = IP.read_depth, cellularity = IP.cellularity)

    return Simulation(IP,simresults,sampleddata)
end


"""
    branchingprocess(IP::InputParameters, rng::AbstractRNG; <keyword arguments>)

Simulate a stochastic branching process with parameters defined by IP. Simulation is by a 
rejection-kinetic Monte Carlo algorithm. If suppressmut assume there is only a single 
acquired mutation per cell at division and no clonal mutations (additional mutations can be 
added retrospectively).

Start simulation with a single cell.

"""
function branchingprocess(IP::InputParameters{BranchingInput}, rng::AbstractRNG,
        fixedmu = IP.fixedmu, μ = IP.μ, clonalmutations = IP.clonalmutations)

    return branchingprocess(IP.siminput.b, IP.siminput.d, IP.siminput.Nmax, μ, rng, 
                            numclones = IP.siminput.numclones, fixedmu = fixedmu, 
                            clonalmutations = clonalmutations, 
                            selection = IP.siminput.selection, tevent = IP.siminput.tevent, 
                            maxclonesize = IP.siminput.maxclonesize, 
                            timefunction = IP.siminput.timefunction)
end

function branchingprocess(b, d, Nmax, μ, rng::AbstractRNG; numclones = 0, fixedmu = false,
    clonalmutations = μ, selection = Float64[], tevent = Float64[], maxclonesize = 200, 
    timefunction::Function = exptime)

    #initialize arrays and parameters
    simtracker = initializesim_branching(b, d, Nmax, rng, numclones = numclones, 
        selection = selection, clonalmutations = clonalmutations)
    
    #run simulation
    simtracker = branchingprocess!(simtracker, Nmax, μ, rng, numclones = numclones,
        fixedmu = fixedmu, tevent = tevent, timefunction = timefunction, 
        maxclonesize = maxclonesize)
    return simtracker
end

"""
    branchingprocess!(simtracker::SimulationTracker, IP::InputParameters, rng::AbstractRNG; 
        suppressmut::Bool = true)

Run branching process simulation, starting in state defined by simtracker, with parameters
defined by IP.

"""
function branchingprocess!(simtracker::BranchingTracker, Nmax, μ, rng::AbstractRNG; 
    numclones = 0, fixedmu = false, tevent = Float64[], maxclonesize = 200, 
    timefunction::Function = exptime)

    t, N = simtracker.tvec[end], simtracker.Nvec[end]
    mutID = N == 1 ? 1 : getmutID(simtracker.cells)

    nclonescurrent = length(simtracker.clonetime) + 1  
    executed = false
    changemutrate = BitArray(undef, numclones + 1)
    changemutrate .= 1

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures
    #that b and d have correct units
    Rmax = (maximum(simtracker.birthrates[1:nclonescurrent])
                                + maximum(simtracker.deathrates[1:nclonescurrent]))

    while N < Nmax
        Nt = N
        randcell = rand(rng,1:N) #pick a random cell
        r = rand(rng,Uniform(0,Rmax))
        #get birth and death rates for randcell
        br = simtracker.birthrates[simtracker.cells[randcell].fitness]
        dr = simtracker.deathrates[simtracker.cells[randcell].fitness]

        if r < br 
            #cell divides
            simtracker, N, mutID = celldivision!(simtracker, randcell, N, mutID, μ, rng,
                fixedmu = fixedmu)
            #check if t>=tevent for next fit subclone
            if nclonescurrent < numclones +1 && t >= tevent[nclonescurrent]
                #if current number clones != final number clones, one of the new cells is
                #mutated to become fitter and form a new clone
                    simtracker, nclonescurrent = cellmutation!(simtracker, N, N, t, 
                        nclonescurrent)
                
                    #change Rmax now there is a new fitter mutant
                    Rmax = (maximum(simtracker.birthrates[1:nclonescurrent])
                                + maximum(simtracker.deathrates[1:nclonescurrent]))
            end

        elseif r < br + dr
            #cell dies
            simtracker, N = celldeath!(simtracker,randcell,N)
        end

        #add time
        Δt =  1/(Rmax * Nt) .* timefunction(rng)
        t = t + Δt
        push!(simtracker.tvec,t)
        #add new population size
        push!(simtracker.Nvec, N)

        if (executed == false) && ((simtracker.clonesize.>maxclonesize) == changemutrate)
            #if population of all clones is sufficiently large no new mutations
            #are acquired, can use this approximation as only mutations above 1%
            #frequency can be reliably detected
            μ = 0
            fixedmu = true
            executed = true
        end
    end
    return simtracker
end

function moranprocess(IP::InputParameters{MoranInput}, rng::AbstractRNG,
                    fixedmu = IP.fixedmu, μ = Ip.μ, clonalmutations = IP.clonalmutations)

    return moranprocess(IP.siminput.N, IP.siminput.bdrate, IP.siminput.tmax, μ, rng, 
                        numclones = IP.siminput.numclones, fixedmu = fixedmu, 
                        clonalmutations = clonalmutations, 
                        selection = IP.siminput.selection, tevent = IP.siminput.tevent, 
                        timefunction = IP.siminput.timefunction)
end

function moranprocess(N, bdrate, tmax, μ, rng::AbstractRNG; numclones = 0, fixedmu = false,
    clonalmutations = μ, selection = Float64[], tevent = Float64[], 
    timefunction::Function = exptime)

    simtracker = initializesim_moran(N, rng, numclones = numclones, clonalmutations = clonalmutations)

    #run simulation
    simtracker = moranprocess!(simtracker, bdrate, tmax, μ, rng, numclones = numclones,
    fixedmu = fixedmu, selection = selection, tevent = tevent, timefunction = timefunction)
    return simtracker
end
"""
    moranprocess!(simtracker::SimulationTracker, rng::AbstractRNG; 
        suppressmut::Bool = true)

Run branching process simulation, starting in state defined by simtracker, with parameters
defined by IP.

"""
function moranprocess!(simtracker::MoranTracker, bdrate, tmax, μ, rng::AbstractRNG; 
    numclones = 0, fixedmu = false, selection = Float64[], tevent = Float64[], 
    timefunction::Function = exptime)

    t, N = simtracker.tvec[end], simtracker.N
    mutID = getmutID(simtracker.cells)

    nclonescurrent = length(simtracker.clonetime) + 1  

    while t < tmax
        #pick a cell to divide proportional to fitness
        if nclonescurrent == 1
            randcelldivide = rand(rng, 1:N)
        else
            p = [cell.fitness==1 ? 1 : 1 + selection[cell.fitness - 1] 
                    for cell in simtracker.cells] 
            p /= sum(p)
            randcelldivide = sample(rng, 1:N, ProbabilityWeights(p)) 
        end

        #pick a random cell to die
        randcelldie = rand(rng,1:N) 

        #cell divides
        simtracker, _, mutID = celldivision!(simtracker, randcelldivide, N, mutID, μ, rng,
            fixedmu = fixedmu)
        
        #check if t>=tevent for next fit subclone and more subclones are expected
        if nclonescurrent < numclones + 1 && t >= tevent[nclonescurrent]
            simtracker, nclonescurrent = cellmutation!(simtracker, N+1, N, t, nclonescurrent)
        end

        #cell dies
        simtracker, _ = celldeath!(simtracker,randcelldie,N)

        #add time
        Δt =  1/(bdrate*N) .* timefunction(rng)
        t = t + Δt
        push!(simtracker.tvec,t)

    end
    return simtracker
end

function set_subclone_birthdeath_rates(b, d, selection, numclones, rng::AbstractRNG)

    birthrates = [b]
    deathrates = [d]
    #add birth and death rates for each subclone. 
    #fitness is randomly distributed between death and birth rates
    for i in 1:numclones
        push!(deathrates, rand(rng) * deathrates[1])
        push!(birthrates,(1 + selection[i]) * (birthrates[1] - deathrates[1]) + deathrates[i + 1])
    end
    return birthrates,deathrates
end

"""
    initializesim(IP::InputParameters, rng::AbstractRNG = Random.GLOBAL_RNG;
    suppressmut = false)

Set up the variables used to track a simulation. If suppressmut is true we assign 0 clonal 
mutations, rather than taking IP.clonalmutations (mutations can be added retrospectively). 

"""
function initializesim_branching(b, d, Nmax, rng::AbstractRNG; numclones = 0, clonalmutations = μ, 
        selection = Float64[])

    birthrates, deathrates = set_subclone_birthdeath_rates(b, d, selection, numclones, rng)

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #population starts with one cell
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population with selection=0
    cells = Cell[]
    sizehint!(cells, Nmax)
    push!(cells, Cell([],1))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    cells[1],mutID = newmutations!(cells[1], clonalmutations, mutID, rng, fixedmu = true)

    clonetype = Int64[] #parent type of each fit subclone 
    clonetime = Float64[]
    acquiredmutations = Array{Int64,1}[]
    cloneN = Int64[]
    Ndivisions = Int64[]
    avdivisions = Float64[]

    clonesize = zeros(Int64, numclones + 1)
    clonesize[1] = 1

    simtracker = BranchingTracker(
        Nvec,
        tvec,
        cells,
        birthrates,
        deathrates,
        clonesize,
        clonetype,
        clonetime,
        acquiredmutations,
        cloneN,
        Ndivisions,
        avdivisions
    )
    return simtracker
end

function initializesim_moran(N, rng::AbstractRNG; numclones = 0, clonalmutations = μ, 
    selection = Float64[])

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population with selection=0
    cells = [Cell([],1) for i in 1:N]

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    for cell in cells
        cell,mutID = newmutations!(cell, clonalmutations, mutID, rng, fixedmu = true)
    end

    clonetype = Int64[] #parent type of each fit subclone 
    clonetime = Float64[]
    acquiredmutations = Array{Int64,1}[]
    cloneN = Int64[]
    Ndivisions = Int64[]
    avdivisions = Float64[]

    clonesize = zeros(Int64, numclones + 1)
    clonesize[1] = N

    simtracker = MoranTracker(
        N,
        tvec,
        cells,
        clonesize,
        clonetype,
        clonetime,
        acquiredmutations,
        cloneN,
        Ndivisions,
        avdivisions
    )
    return simtracker
end

function newmutations!(cell, μ, mutID, rng::AbstractRNG; fixedmu = true)
    #function to add new mutations to cells based on μ
    numbermutations = fixedmu ? μ : rand(rng,Poisson(μ))
    append!(cell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations
    return cell, mutID
end

function celldivision!(simtracker::SimulationTracker, parentcell, N, mutID, μ, 
    rng::AbstractRNG; fixedmu=fixedmu)
    
    N += 1 #population increases by one
    push!(simtracker.cells, copycell(simtracker.cells[parentcell])) #add new copy of parent cell to cells
    #add new mutations to both new cells
    if μ > 0.0 
        simtracker.cells[parentcell],mutID = newmutations!(simtracker.cells[parentcell], μ, 
            mutID, rng, fixedmu = fixedmu)
        simtracker.cells[end],mutID = newmutations!(simtracker.cells[end], μ, mutID, rng,
            fixedmu=fixedmu)
    end
    simtracker.clonesize[simtracker.cells[parentcell].fitness] += 1
    return simtracker, N, mutID
end

function cellmutation!(simtracker::SimulationTracker, mutatingcell, N, t, nclonescurrent)

    simtracker.clonesize[simtracker.cells[mutatingcell].fitness] -= 1
    nclonescurrent += 1
    push!(simtracker.clonetype, simtracker.cells[mutatingcell].fitness) 

    #change fitness type of new cell
    simtracker.cells[mutatingcell].fitness = nclonescurrent
    simtracker.clonesize[simtracker.cells[mutatingcell].fitness] += 1

    push!(simtracker.clonetime, t)
    push!(simtracker.acquiredmutations, deepcopy(simtracker.cells[mutatingcell].mutations))
    push!(simtracker.cloneN, N)
    push!(simtracker.Ndivisions, length(simtracker.cells[mutatingcell].mutations))
    divs = map(x -> length(x.mutations), simtracker.cells)
    push!(simtracker.avdivisions, mean(divs))

    return simtracker, nclonescurrent
end

function celldeath!(simtracker::SimulationTracker, deadcell, N)
    #population decreases by 1
    N -= 1
    #frequency of cell type decreases
    simtracker.clonesize[simtracker.cells[deadcell].fitness] -= 1
    #remove deleted cell
    deleteat!(simtracker.cells,deadcell)
    return simtracker, N
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
    return Cell(copy(cellold.mutations), copy(cellold.fitness))
  end

function exptime(rng::AbstractRNG)
    - log(rand(rng))
end

function no_mutations(cell)
    return length(cell.mutations) == 0
end
