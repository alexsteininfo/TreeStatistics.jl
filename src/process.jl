function processresults!(moduletracker::ModuleTracker, μ, clonalmutations, rng::AbstractRNG;
    mutationdist=:poisson)

    mutationlist = get_mutationlist(moduletracker)
    expandedmutationids = 
        get_expandedmutationids(μ, mutationlist, clonalmutations, rng, mutationdist=mutationdist)
    expandmutations!(moduletracker, expandedmutationids, clonalmutations)
    
    return moduletracker
end

function processresults!(populationtracker::Array{ModuleTracker, 1}, μ, clonalmutations, 
    rng::AbstractRNG; mutationdist=:poisson)
    
    mutationlist = get_mutationlist(populationtracker)
    expandedmutationids = 
        get_expandedmutationids(μ, mutationlist, clonalmutations, rng, mutationdist=mutationdist)
    
    for moduletracker in populationtracker
        expandmutations!(moduletracker, expandedmutationids, clonalmutations)
    end
    return populationtracker
end

function final_timedep_mutations!(populationtracker::Array{ModuleTracker, 1}, μ, mutationdist, rng)
    mutID = maximum(mutid 
        for moduletracker in populationtracker 
            for cell in moduletracker.cells
                for mutid in cell.mutations
    )
    tend = age(populationtracker)
    for moduletracker in populationtracker
        mutID = final_timedep_mutations!(moduletracker, μ, mutationdist, tend, rng, mutID=mutID)
    end
end

function final_timedep_mutations!(moduletracker::ModuleTracker, μ, mutationdist, rng; mutID=nothing)
    tend = age(moduletracker)
    if isnothing(mutID)
        mutID = maximum(mutid for cell in moduletracker.cells for mutid in cell.mutations) + 1
    end
    for cell in moduletracker.cells
        Δt = tend - cell.birthtime
        numbermutations = numbernewmutations(rng, mutationdist, μ, Δt=Δt)
        mutID = addnewmutations!(cell, numbermutations, mutID)
    end
    return mutID
end

function expandmutations!(moduletracker, expandedmutationids, clonalmutations)

    if length(expandedmutationids) > 0
        for cell in moduletracker.cells
            cell.mutations = expandmutations(expandedmutationids, cell.mutations)
            if clonalmutations > 0
                prepend!(cell.mutations, 1:clonalmutations)
            end
        end
    elseif clonalmutations > 0
        for cell in moduletracker.cells
            cell.mutations = collect(1:clonalmutations)
        end
    end
    #get list of mutations in each subclone
    if length(expandedmutationids) > 0
        for subclone in moduletracker.subclones
            subclone.mutations = expandmutations(expandedmutationids, subclone.mutations)
        end
    end
    return moduletracker
end

function expandmutations(expandedmutationids, originalmutations)
    return reduce(
        vcat, 
        filter(!isempty, map(x -> expandedmutationids[x], originalmutations)),
        init=Int64[]
    )
end

function get_mutationlist(populationtracker::Array{ModuleTracker, 1})
    #get list of all mutations assigned to each cell
    mutationlist = [mutation 
        for moduletracker in populationtracker
            for cell in moduletracker.cells
                for mutation in cell.mutations
    ]
    return sort(unique(mutationlist))
end

function get_mutationlist(moduletracker::ModuleTracker)
    #get list of all mutations assigned to each cell
    mutationlist = [mutation 
        for cell in moduletracker.cells
            for mutation in cell.mutations
    ]
    return unique(mutationlist)
end

function get_expandedmutationids(μ, mutationlist, clonalmutations, rng; mutationdist=:poisson)

    mutationsN = numberuniquemutations(rng, length(mutationlist), mutationdist, μ)
    expandedmutationids = Dict{Int64, Vector{Int64}}()
    i = clonalmutations + 1
    for (mutkey, N) in zip(mutationlist, mutationsN)
        push!(expandedmutationids, mutkey=>collect(i:i+N-1))
        i += N
    end
    return expandedmutationids
end

function numberuniquemutations(rng, L, mutationdist, μ)
    #returns the number of unique mutations drawn from a distribution according to the
    #mutation rule with mean μ, L times.
    if mutationdist == :fixed
        return fill(μ, L)
    elseif mutationdist == :poisson
        return rand(rng, Poisson(μ), L) 
    elseif mutationdist == :geometric
        return rand(rng, Geometric(1/(1+μ)), L)
    else
        error("$mutationdist is not a valid mutation rule")
    end
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
