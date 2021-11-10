mutable struct Cell
    mutations::Array{Int64,1}
    fitness::Int64
end

struct SimulationTracker
    Nvec::Array{Int64, 1}
    tvec::Array{Float64, 1}
    muts::Array{Int64, 1}
    cells::Array{Cell, 1}
    birthrates::Array{Float64, 1}
    deathrates::Array{Float64, 1}
    clonesize::Array{Int64,1}
    clonetype::Array{Int64, 1}
    clonetime::Array{Float64, 1}
    acquiredmutations::Array{Array{Int64, 1}, 1}
    cloneN::Array{Int64, 1}
    Ndivisions::Array{Int64, 1}
    avdivisions::Array{Float64, 1}
end

abstract type SimulationInput end

struct BranchingInput <: SimulationInput
    numclones::Int64 
    Nmax::Int64
    clonalmutations::Int64
    selection::Array{Float64,1}
    μ::Float64
    b::Float64
    d::Float64
    tevent::Array{Float64,1}
    fixedmu::Bool
    timefunction::Function
    maxclonesize::Int64
end

struct MoranInput <: SimulationInput
    N::Int64
    numclones::Int64 
    tmax::Float64
    clonalmutations::Int64
    selection::Array{Float64,1}
    μ::Float64
    bdrate::Float64
    tevent::Array{Float64,1}
    fixedmu::Bool
    timefunction::Function
    maxclonesize::Int64
end

struct InputParameters{T<:SimulationInput}
    detectionlimit::Float64
    ploidy::Int64
    read_depth::Float64
    ρ::Float64
    cellularity::Float64
    siminput::T
end

struct SimulationResult
    clonefreq::Array{Float64,1}
    clonefreqp::Array{Float64,1}
    clonetime::Array{Float64,1}
    subclonalmutations::Array{Int64,1}
    birthrates::Array{Float64,1}
    deathrates::Array{Float64,1}
    tend::Float64
    trueVAF::Array{Float64,1}
    cloneN::Array{Int64, 1}
    clonetype::Array{Int64, 1}
    Ndivisions::Array{Int64, 1}
    cells::Array{Cell, 1}
    avdivisions::Array{Float64, 1}
end

struct SampledData
    df::DataFrame
    VAF::Array{Float64,1}
    counts::Array{Int64,1}
    depth::Array{Int64,1}
end

struct Simulation{T<:SimulationInput}
    input::InputParameters{T}
    output::SimulationResult
    sampled::SampledData
end

function InputParameters{BranchingInput}(;numclones = 1, Nmax = 10000, ploidy = 2, 
    read_depth = 100.0, detectionlimit = 5/read_depth, μ = 10.0, clonalmutations = μ, 
    selection = fill(0.0,numclones), b = log(2.0), d = 0.0, 
    tevent = collect(1.0:0.5:(1+numclones)/2), ρ = 0.0, cellularity = 1.0, fixedmu = false, 
    timefunction = exptime, maxclonesize = 200)

    return InputParameters(
        detectionlimit,
        ploidy,
        read_depth,
        ρ,
        cellularity,
        BranchingInput(
            numclones,
            Nmax,
            clonalmutations,
            selection,
            μ,
            b,
            d,
            tevent,
            fixedmu,
            timefunction,
            maxclonesize
        )
    )
end

# function InputParametersMoran(IP::InputParameters,tmax,bdrate)

#     return InputParametersMoran(
#         IP.numclones,
#         tmax,
#         IP.detectionlimit,
#         IP.ploidy,
#         IP.read_depth,
#         IP.clonalmutations,
#         IP.selection,
#         IP.μ,
#         bdrate,
#         IP.tevent,
#         IP.ρ,
#         IP.cellularity,
#         IP.fixedmu,
#         IP.timefunction,
#         IP.maxclonesize
#     )
# end

"""
    run1simulation(b, d, Nmax; <keyword arguments>)
"""
function run1simulation(IP::InputParameters{BranchingInput}, rng::AbstractRNG = MersenneTwister();
    minclonefreq = 0.0, maxclonefreq = 1.0)

    #Initially set clonalmutations = 0 and μ=1. These are expanded later.
    simtracker = 
        branchingprocess(IP.siminput.b, IP.siminput.d, IP.siminput.Nmax, 1, rng, numclones = IP.siminput.numclones, 
            fixedmu = true, clonalmutations = 0, selection = IP.siminput.selection,
            tevent = IP.siminput.tevent, maxclonesize = IP.siminput.maxclonesize, 
            timefunction = IP.siminput.timefunction)
        
    simtracker, numclones, simresults = 
        processresults!(simtracker, IP.siminput.Nmax, IP.siminput.numclones, IP.siminput.μ, 
                        IP.siminput.fixedmu, IP.siminput.clonalmutations, IP.ploidy, 
                        minclonefreq, maxclonefreq, rng)
    
    sampleddata = 
        sampledhist(simresults.trueVAF, IP.siminput.Nmax, rng, 
                    detectionlimit = IP.detectionlimit, 
                    read_depth = IP.read_depth, cellularity = IP.cellularity)

    return Simulation(IP,simresults,sampleddata)
end

function run1simulation(IP::InputParameters{MoranInput}, rng::AbstractRNG = MersenneTwister();
    minclonefreq = 0.0, maxclonefreq = 1.0)

    #Initially set clonalmutations = 0 and μ=1. These are expanded later.
    simtracker = 
        moranprocess(IP.siminput.b, IP.siminput.d, IP.siminput.Nmax, 1, rng, 
                    numclones = IP.siminput.numclones, fixedmu = true, 
                    clonalmutations = 0, selection = IP.siminput.selection,
                    tevent = IP.siminput.tevent, maxclonesize = IP.siminput.maxclonesize, 
                    timefunction = IP.siminput.timefunction)
        
    simtracker, numclones, simresults = 
        processresults!(simtracker, IP.siminput.Nmax, IP.siminput.numclones, IP.siminput.μ, 
                        IP.siminput.fixedmu, IP.siminput.clonalmutations, 
                        IP.ploidy, minclonefreq, maxclonefreq, rng)
    
    sampleddata = 
        sampledhist(simresults.trueVAF, IP.siminput.Nmax, rng, 
                    detectionlimit = IP.detectionlimit, 
                    read_depth = IP.read_depth, cellularity = IP.cellularity)

    return SimulationBranching(IP,simresults,sampleddata)
end


"""
    branchingprocess(IP::InputParameters, rng::AbstractRNG; suppressmut::Bool = true)

Simulate a stochastic branching process with parameters defined by IP. Simulation is by a 
rejection-kinetic Monte Carlo algorithm. If suppressmut assume there is only a single 
acquired mutation per cell at division and no clonal mutations (additional mutations can be 
added retrospectively).

Start simulation with a single cell.

"""
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
function branchingprocess!(simtracker::SimulationTracker, Nmax, μ, rng::AbstractRNG; 
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
        bcell = simtracker.birthrates[simtracker.cells[randcell].fitness]
        dcell = simtracker.deathrates[simtracker.cells[randcell].fitness]

        if r < bcell 
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

        elseif r < bcell + dcell
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

function moranprocess(N, bdrate, tmax, μ, rng::AbstractRNG; numclones = 0, fixedmu = false,
    clonalmutations = μ, selection = Float64[], tevent = Float64[], maxclonesize = 200, 
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
function moranprocess!(simtracker::SimulationTracker, bdrate, tmax, μ, rng::AbstractRNG; 
    numclones = 0, fixedmu = false, selection = Float64[], tevent = Float64[], 
    timefunction::Function = exptime)

    t, N = simtracker.tvec[end], simtracker.Nvec[end]
    mutID = getmutID(simtracker.cells)

    nclonescurrent = length(simtracker.clonetime) + 1  

    while t < tmax
        N0 = N
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
        simtracker, N, mutID = celldivision!(simtracker, randcelldivide, N, mutID, μ, rng,
            fixedmu = fixedmu)
        
        #check if t>=tevent for next fit subclone and more subclones are expected
        if nclonescurrent < numclones + 1 && t >= tevent[nclonescurrent]
            simtracker, nclonescurrent = cellmutation!(simtracker, N, N0, t, nclonescurrent)
        end

        #cell dies
        simtracker, N = celldeath!(simtracker,randcelldie,N)

        #add time
        Δt =  1/bdrate .* timefunction(rng)
        t = t + Δt
        push!(simtracker.tvec,t)
        #add new population size
        push!(simtracker.Nvec, N)

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
    initializesim(IP::InputParameters, rng::AbstractRNG = MersenneTwister();
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
    muts = Int64[]
    push!(muts,mutID)

    clonetype = Int64[] #parent type of each fit subclone 
    clonetime = Float64[]
    acquiredmutations = Array{Int64,1}[]
    cloneN = Int64[]
    Ndivisions = Int64[]
    avdivisions = Float64[]

    clonesize = zeros(Int64, numclones + 1)
    clonesize[1] = 1

    simtracker = SimulationTracker(
        Nvec,
        tvec,
        muts,
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

    #population starts with one cell
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population with selection=0
    cells = [Cell([],1) for i in 1:N]

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    for cell in cells
        cell,mutID = newmutations!(cell, clonalmutations, mutID, rng, fixedmu = true)
        muts = Int64[]
        push!(muts,mutID)
    end

    clonetype = Int64[] #parent type of each fit subclone 
    clonetime = Float64[]
    acquiredmutations = Array{Int64,1}[]
    cloneN = Int64[]
    Ndivisions = Int64[]
    avdivisions = Float64[]

    clonesize = zeros(Int64, numclones + 1)
    clonesize[1] = 1

    simtracker = SimulationTracker(
        Nvec,
        tvec,
        muts,
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
    push!(simtracker.muts,mutID)
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
    allmutations = reduce(vcat,[cell.mutations for cell in cells])
    return maximum(allmutations)
end

function copycell(cellold::Cell)
    return Cell(copy(cellold.mutations), copy(cellold.fitness))
  end

function exptime(rng::AbstractRNG)
    - log(rand(rng))
end

# function meantime()
#     1
# end
