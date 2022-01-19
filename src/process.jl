function processresults!(moduletracker::ModuleTracker, μ, clonalmutations, rng::AbstractRNG)

    maxmutation = maximum(maximum(cell.mutations) for cell in moduletracker.cells)
    mutationsN = rand(rng, Poisson(μ), maxmutation) 
    mutationids = get_mutationids(μ, maxmutation, clonalmutations, rng)
    
    for cell in moduletracker.cells
        cell.mutations = reduce(vcat, mutationids[cell.mutations])
        if clonalmutations > 0
            prepend!(cell.mutations, 1:clonalmutations)
        end
    end

    #get list of mutations in each subclone
    for subclone in moduletracker.subclones
        subclone.mutations = reduce(vcat, mutationids[subclone.mutations])
        prepend!(subclone.mutations, 1:clonalmutations)
    end

    return moduletracker
end

function processresults!(populationtracker::Array{ModuleTracker, 1}, μ, clonalmutations, 
    rng::AbstractRNG)
    
    maxmutation = maximum(maximum(cell.mutations) 
        for moduletracker in populationtracker for cell in moduletracker.cells
    )
    mutationids = get_mutationids(μ, maxmutation, clonalmutations, rng)
    
    for moduletracker in populationtracker
        for cell in moduletracker.cells
            cell.mutations = reduce(vcat, mutationids[cell.mutations])
            if clonalmutations > 0
                prepend!(cell.mutations, 1:clonalmutations)
            end
        end

        #get list of mutations in each subclone
        for subclone in moduletracker.subclones
            subclone.mutations = reduce(vcat, mutationids[subclone.mutations])
            prepend!(subclone.mutations, 1:clonalmutations)
        end
    end
    return populationtracker
end

function get_mutationids(μ, maxmutation, clonalmutations, rng)
    mutationsN = rand(rng, Poisson(μ), maxmutation) 
    mutationids = Vector{Int64}[]
    i = clonalmutations + 1
    for N in mutationsN
        push!(mutationids, collect(i:i+N-1))
        i += N
    end
    return mutationids
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
