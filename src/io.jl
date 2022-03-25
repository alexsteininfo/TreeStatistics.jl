typedict(x) = Dict(fn=>getfield(x, fn) for fn in fieldnames(typeof(x)))

function saveinput(input, filename)
    filename = filename[end-4:end] == ".json" ? filename : filename * ".json"
    inputdict = typedict(input)
    open(filename, "w") do io
        JSON.print(io, inputdict, 4)
    end
end

function loadinput(inputtype, filename)
    inputdict = JSON.parsefile(filename, dicttype=Dict{Symbol, Any})
    loadinput(inputtype, inputdict::Dict)
end

function loadinput(inputtype, inputdict::Dict)
    return inputtype(;inputdict...)
end


function saveoutput(population, output, outputdir, seed, id, rng::AbstractRNG=Random.GLOBAL_RNG)
    outputdir = outputdir[end] == '/' ? outputdir : outputdir * '/'
    mkpath(outputdir)

    if haskeytrue(output, :compare_sampled_modules)
        save_sampled_module_compare(population, output[:nsample], outputdir, id, rng)
    end
    if haskeytrue(output, :pfd_matrix)
        save_pfd_matrix(population, outputdir, id)
    end
    if haskeytrue(output, :sharedfixedmuts)
        save_sharedfixedmuts(population, outputdir, id)
    end
    if haskeytrue(output, :moduleancestory)
        save_moduleancestory(population, outputdir, id)
    end
    save_other_population_data(population, output, outputdir, id)

end

haskeytrue(output, key) = haskey(output, key) && output[key]    

function save_moduleancestory(population, outputdir, id)
    open(outputdir*"moduleancestory.txt", "w") do io
        for moduletracker in population
            write(io, "$(moduletracker.parentid)    $(moduletracker.id)    $(moduletracker.tvec[1])\n")
        end
    end
end

function save_sampled_module_compare(population, nsample, outputdir, id, rng)
    pairwisefixeddiff, sharedfixedmuts = sampledmoduledata(population, nsample, rng)
        open(outputdir*"pairwisefixeddiff_$id.txt", "w") do io
            writedlm(io, zip(keys(pairwisefixeddiff), values(pairwisefixeddiff)))
        end
        open(outputdir*"sharedfixedmuts_$id.txt", "w") do io
            writedlm(io, zip(keys(sharedfixedmuts), values(sharedfixedmuts)))
        end
end

function sampledmoduledata(population, nsample, rng)
    samplemodules = sample(rng, 1:length(population), nsample, replace=false)
    clonalmuts = clonal_mutation_ids(population, samplemodules)
    pairwisefixeddiff = pairwise_fixed_differences(clonalmuts)
    sharedfixedmuts = shared_fixed_mutations(clonalmuts)
    return pairwisefixeddiff, sharedfixedmuts
end

function save_pfd_matrix(population, outputdir, id)
    pfd = pairwise_fixed_differences_matrix(population, diagonals=true)
    open(outputdir*"pfdmatrix_$id.txt", "w") do io
        writedlm(io, pfd)
    end
end

function save_sharedfixedmuts(population, outputdir, id)
    sharedfixedmuts = shared_fixed_mutations(population)
    open(outputdir*"sharedfixedmutsall_$i.txt", "w") do io
        writedlm(io, sharedfixedmuts)
    end
end

function save_other_population_data(population, output, outputdir, id)
    outputdict = Dict()
    if output[:finalmodules]
        outputdict[:finalmodules] = length(population)
    end
    if output[:finalage]
        outputdict[:finalage] = age(population)
    end
    if output[:mutation_stats]
        outputdict[:mutations_mean], outputdict[:mutations_var] = 
            average_mutations(population, true)
    end
    if output[:pfd_stats]
        pfd_mean, pfd_var, fixedmutations_mean, fixedmutations_var = 
            pairwise_fixed_differences_statistics(population, clonal=true)
        outputdict[:pfd_mean] = pfd_mean
        outputdict[:pfd_var] = pfd_var
        outputdict[:fixedmutations_mean] = fixed_mutations_mean
        outputdict[:fixedmutations_var] = fixed_mutations_var
    end
    open(outputdir*"summarydata_$id.json", "w") do io
        JSON.print(io, outputdict, 4)
    end
end