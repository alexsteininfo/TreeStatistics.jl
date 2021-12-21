function processresults!(moduletracker, Nmax, numclones, μ, fixedmu::Bool, 
    clonalmutations, ploidy, rng::AbstractRNG)

    if length(moduletracker.subclones)!= numclones
        error("wrong number of clones, $(length(moduletracker.subclones)) != $numclones")
    end

    #get list of all mutations (assumed one per cell division)
    mutations = cellsconvert(moduletracker.cells).mutations
    #get list of mutations in each subclone
    subclonalmutation_ids = [subclone.mutations for subclone in moduletracker.subclones]

    #calculate allele frequencies. Expand so that each cell obtained m ~ Poisson(μ) 
    #mutations at division and add clonal mutations
    allelefreq = getallelefreq(mutations, Nmax)
    allelefreq, subclonalmutations = 
        allelefreqexpand(allelefreq, μ, subclonalmutation_ids, rng, fixedmu = fixedmu)
    prepend!(allelefreq, fill(Float64(Nmax), clonalmutations))
    allelefreq /= (ploidy * Nmax) #correct for ploidy
    
    #get proportion of cells in each subclone
    clonesize = getclonesize(moduletracker)
    clonefreqp = clonesize[2:end]/sum(clonesize)
    clonefreq = copy(clonefreqp)
    if length(clonefreqp) > 1
        clonefreq, subclonalmutations = 
            calculateclonefreq!(clonefreq, subclonalmutations, moduletracker.subclones)
    end

    subclones = [
        Clone(
            subclone.parenttype, 
            subclone.parentmodule, 
            subclone.time, 
            submuts, 
            subclone.N0, 
            subclone.Ndivisions, 
            subclone.avdivisions, 
            freq, 
            freqp
        )
        for (subclone, submuts, freq, freqp) 
            in zip(moduletracker.subclones, subclonalmutations, clonefreq, clonefreqp)]
    

    return moduletracker, SimulationResult(
        subclones,
        moduletracker.tvec[end], 
        allelefreq, 
        moduletracker.cells,
        moduletracker.Nvec,
        moduletracker.tvec
        )   
end


function get_pophistory(moduletrackervec::Array{ModuleTracker, 1})
    Nvec, tvec = get_pophistory(moduletrackervec[1])
    for moduletracker in moduletrackervec[2:end]
        Nvec0, tvec0 = get_pophistory(moduletracker)
        append!(Nvec, Nvec0[2:end])
        append!(tvec, tvec0[2:end])
    end
    return Nvec, tvec
end

function cellsconvert(cells)
    #convert from array of cell types to one array with mutations and one array with cell fitness

    clonetype = zeros(Int64,length(cells))
    mutations = Int64[]
    # sizehint!(mutations, length(cells) * 10) #helps with performance to provide vector size

    for i in 1:length(cells)
        append!(mutations,cells[i].mutations)
        clonetype[i] = cells[i].clonetype
    end

    return (;mutations, clonetype)
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

function allelefreqexpand(AFDict, μ, subclonalmuation_ids, rng::AbstractRNG; fixedmu = false)

    #expand allele frequency given mutation rate and calculate number of mutations in the subclones
    
    subclonalmutations = zeros(Int64, length(subclonalmuation_ids))
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
        idx = findall((in)(subclonalmuation_ids[i]), mutids) 
        subclonalmutations[i] = sum(mutations[idx]) 
    end

    j = 0
    for f in 1:length(mutfreqs)
        AFnew[(j + 1): j + mutations[f]] = fill(mutfreqs[f], mutations[f])
        j = j + mutations[f]
    end
  return AFnew, subclonalmutations
end

function remove_undetectable!(moduletracker::ModuleTracker, clonefreq, clonefreqp, numclones, detectableclones)
    #if there are clones outside the detectable range remove them from the data
    if sum(detectableclones) < numclones
        numclones = sum(detectableclones)
        clonefreq = clonefreq[detectableclones]
        clonefreqp = clonefreqp[detectableclones]
        moduletracker.subclones = moduletracker.subclones[detectableclones]
        pushfirst!(detectableclones, true)
        moduletracker.clonesize = moduletracker.clonesize[detectableclones]
        detectableclones = detectableclones[1:length(br)]
    end
    return moduletracker, clonefreq, clonefreqp, numclones
end


function calculateclonefreq!(clonefreq, subclonalmutations, subclones)
    for i in length(subclones):-1:2
        if subclones[i].parenttype > 1
            clonefreq[subclones[i].parenttype-1] += clonefreq[i]
            subclonalmutations[i] -= subclonalmutations[subclones[i].parenttype-1]
        end
    end
    if (sum(clonefreq.>1.0) > 0)
        error("There is a clone with frequency greater than 1, this should be impossible 
                ($(clonesize)), $(parenttype), $(clonefreq)")
    end
    return clonefreq, subclonalmutations
  end
