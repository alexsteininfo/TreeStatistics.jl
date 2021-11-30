function processresults!(simtracker::SimulationTracker, Nmax, numclones, μ, fixedmu::Bool, 
    clonalmutations, ploidy, minclonefreq, maxclonefreq, rng::AbstractRNG)

    mutations = cellsconvert(simtracker.cells).mutations
    
    if length(simtracker.clonetime)!= numclones
        error("wrong number of clones")
    end

    allelefreq = getallelefreq(mutations, Nmax)
    allelefreq, subclonalmutations = 
        allelefreqexpand(allelefreq, μ, simtracker.acquiredmutations, rng, 
                        fixedmu = fixedmu)

    prepend!(allelefreq, fill(Float64(Nmax), clonalmutations))
    allelefreq /= (ploidy * Nmax)
    
    #get proportion of cells in each subclone
    clonefreqp = simtracker.clonesize[2:end]/sum(simtracker.clonesize)
    clonefreq = copy(clonefreqp)
    if length(clonefreqp) > 1
        clonefreq, subclonalmutations = 
            calculateclonefreq!(clonefreq, subclonalmutations, simtracker.clonetype)
    end

    detectableclones = (clonefreq .> minclonefreq) .& (clonefreq .< maxclonefreq)
    
    if sum(detectableclones) > numclones
        numclones, clonefreq, clonefreqp, simtracker =
            remove_undetectable!(numclones, clonefreq, clonefreqp, simtracker, 
                                detectableclones)
    end
    
    return simtracker, SimulationResult(
        clonefreq, 
        clonefreqp, 
        simtracker.clonetime,
        subclonalmutations, 
        simtracker.tvec[end], 
        allelefreq, 
        simtracker.cloneN, 
        simtracker.clonetype, 
        simtracker.Ndivisions, 
        simtracker.cells, 
        simtracker.avdivisions)
end


function cellsconvert(cells)
    #convert from array of cell types to one array with mutations and one array with cell fitness

    fitness = zeros(Int64,length(cells))
    mutations = Int64[]
    # sizehint!(mutations, length(cells) * 10) #helps with performance to provide vector size

    for i in 1:length(cells)
        append!(mutations,cells[i].mutations)
        fitness[i] = cells[i].fitness
    end

    return (;mutations, fitness)
end

function getallelefreq(mutations,cellnum)
    #create dictionary that maps mutation ID to allele frequency
    f = counts(mutations,minimum(mutations):maximum(mutations))
    muts = collect(minimum(mutations):maximum(mutations))
    idx = f .> 0.01 #should this be f/(2*cellnum) .> 0.01, i.e. only include freq > 1% ??
    f = map(Float64, f[idx])
    muts = muts[idx]
    Dict{Int64, Float64}(muts[i]::Int64 => f[i]::Float64 for i in 1:length(f))
end

function allelefreqexpand(AFDict, μ, acquiredmutations, rng::AbstractRNG; fixedmu = false)

    #expand allele frequency given mutation rate and calculate number of mutations in the subclones
    #acquiredmutations = convert(Array{Array{Int64,1},1}, acquiredmutations)
    
    subclonalmutations = zeros(Int64, length(acquiredmutations))
    mutfreqs = collect(values(AFDict))
    mutids = collect(keys(AFDict))
    if fixedmu
        μ = round(Int64, μ) #ensure number of mutations is a whole number
        mutations = fill(μ, length(mutfreqs))
    else 
        #randomly sample number of mutations
        mutations = rand(rng, Poisson(μ), length(mutfreqs)) 
    end
    mutations = rand(rng, Poisson(μ), length(mutfreqs)) #randomly sample number of mutations
    AFnew = zeros(Int64, sum(mutations))
    #find all subclonal mutations, i.e. mutations held by all members of each fit subclone
    for i in 1:length(subclonalmutations)
        idx = findall((in)(acquiredmutations[i]), mutids) 
        subclonalmutations[i] = sum(mutations[idx]) 
    end

    j = 0
    for f in 1:length(mutfreqs)
        AFnew[(j + 1): j + mutations[f]] = fill(mutfreqs[f], mutations[f])
        j = j + mutations[f]
    end
  return AFnew, subclonalmutations
end

function remove_undetectable!(numclones, clonefreq, clonefreqp, simtracker, detectableclones)
    #if there are clones outside the detectable range remove them from the data
    if sum(detectableclones) != numclones
        numclones = sum(detectableclones)
        clonefreq = clonefreq[detectableclones]
        clonefreqp = clonefreqp[detectableclones]
        simtracker.clonetime = simtracker.clonetime[detectableclones]
        simtracker.clonetype = simtracker.clonetype[detectableclones]
        pushfirst!(detectableclones, true)
        detectableclones = detectableclones[1:length(br)]
        simtracker.birthrates = simtracker.birthrates[detectableclones]
        simtracker.deathrates = simtracker.deathrates[detectableclones]
    end
    return numclones, clonefreq, clonefreqp, simtracker
end

function calculateclonefreq!(clonefreq, subclonalmutations, parenttype)
    for i in length(parenttype):-1:2
        if parenttype[i] > 1
            clonefreq[parenttype[i]-1] += clonefreq[i]
            subclonalmutations[i] -= subclonalmutations[parenttype[i]-1]
        end
    end
    if (sum(clonefreq.>1.0) > 0)
        error("There is a clone with frequency greater than 1, this should be impossible 
                ($(clonesize)), $(parenttype), $(clonefreq)")
    end
    return clonefreq, subclonalmutations
  end
