
function multiplesimulations(numsim, IPlist... ; rng::AbstractRNG = Random.GLOBAL_RNG,
    minclonefreq = 0.0, maxclonefreq = 1.0)

    multsimlist = [multiplesimulations(numsim, IP, rng = rng, minclonefreq = minclonefreq, 
                                        maxclonefreq = maxclonefreq) 
                        for IP in IPlist]
    return multsimlist
end

function multiplesimulations(numsim, IP; rng::AbstractRNG = Random.GLOBAL_RNG,
    minclonefreq = 0.0, maxclonefreq = 1.0)
    
    results = MultiSimulation(IP, SimulationResult[], SampledData[])
    for i in 1:numsim
        simdata = run1simulation(IP, rng, minclonefreq = minclonefreq,
            maxclonefreq = maxclonefreq)
        push!(results.output, simdata.output)
        push!(results.sampled, simdata.sampled)
    end
    return results
end

function getstats(multsimlist::Array{MultiSimulation{T},1}, args::Symbol...) where T <: SimulationInput
    i = 1
    dflist = DataFrame[]
    for multsim in multsimlist
        μ_predicted = Float64[]
        r2 = Float64[]
        for simresults in multsim.sampled
            _, fitcoef, r2val = fitinverse(simresults.VAF, 0.12, 0.24)
            push!(μ_predicted, fitcoef)
            push!(r2, r2val)
        end 
        μ_error = ((μ_predicted .- multsim.input.siminput.μ) 
                        ./ multsim.input.siminput.μ .* 100)
        df = DataFrame(id = i, μ_predicted = μ_predicted, μ_error = μ_error, r2 = r2)
        addinputparams!(df,multsim,args)
        i += 1
        push!(dflist, df)
    end
    df = reduce(vcat, dflist, cols = :union)
end

function addinputparams!(df,multsim,args)
    for fieldname in args
        if fieldname in fieldnames(InputParameters)
            value = getfield(multsim.input, fieldname)
        elseif fieldname in fieldnames(typeof(multsim.input.siminput))
            value = getfield(multsim.input.siminput, fieldname)
        end
        addparam!(df,fieldname,value)
    end
end

function addparam!(df,fieldname,value)
    if value isa AbstractArray
        for (i,v) in enumerate(value)
            df[!, Symbol(fieldname, i)] .= v
        end
    else
        df[!, fieldname] .= value
    end
    return df
end