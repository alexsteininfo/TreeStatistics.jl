
function multiplesimulations(numsim, inputlist... ; rng::AbstractRNG = Random.GLOBAL_RNG)

    multsimlist = [multiplesimulations(numsim, input, rng = rng) 
                        for input in inputlist]
    return multsimlist
end

function multiplesimulations(::Type{T}, numsim, input; rng::AbstractRNG = Random.GLOBAL_RNG) where T
    
    results = MultiSimulation(input, CellModule[])
    for i in 1:numsim
        simdata = runsimulation(T, input, rng)
        push!(results.output, simdata.output)
    end
    return results
end

function getstats(multsimlist::Vector{MultiSimulation{T}}, args::Symbol...;
    fit_fmin=0.12, fit_fmax=0.24, sampled=true) where T <: SimulationInput
    i = 1
    dflist = DataFrame[]
    for multsim in multsimlist
        μ_predicted = Float64[]
        r2 = Float64[]
        data = sampled ? multsim.sampled : multsim.output
        for simresults in data
            VAF = sampled ? simresults.VAF : simresults.trueVAF
            _, fitcoef, r2val = fitinverse(VAF, fit_fmin, fit_fmax)
            push!(μ_predicted, fitcoef[1])
            push!(r2, r2val)
        end 
        μ_error = ((μ_predicted .- multsim.input.μ) 
                        ./ multsim.input.μ .* 100)
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
        elseif fieldname in fieldnames(typeof(multsim.input))
            value = getfield(multsim.input, fieldname)
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