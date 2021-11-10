
"""
    sampledhist(trueVAF::Array{Float64,1}, ncells::Int64; <keyword arguments>)

    Sample from the true 
"""
function sampledhist(trueVAF::Array{Float64, 1}, cellnum::Int64, rng::AbstractRNG; 
    detectionlimit = 5/read_depth, read_depth = 100.0, cellularity = 1.0)

    VAF = trueVAF * cellularity
    filter!(x -> x > detectionlimit, VAF)
    depth = rand(rng, Binomial(cellnum,read_depth/cellnum), length(VAF))
    # depth = rand(Poisson(read_depth), length(trueVAF)) #why is depth Poisson distributed??
    sampalleles = map((n, p) -> rand(rng, Binomial(n, p)), depth, VAF)

    sampledVAF = sampalleles./depth
    dfhist = gethist(sampledVAF,0.005,0.995,0.01)
    #data for histogram
    return SampledData(dfhist, sampledVAF, sampalleles, depth)
end

function gethist(VAF::Array{Float64,1}, xmin=0.005, xmax=0.995, xstep=0.01)
    x = xmin:xstep:xmax+xstep
    y = fit(Histogram, VAF, x, closed=:right).weights
    dfhist = DataFrame(VAF = x[1:end-1], freq = y)
    dfhist = addcumfreq!(dfhist,:freq)
    return dfhist
end

function addcumfreq!(df,colname)
    df[!, Symbol(:cum,colname)] = reverse(cumsum(reverse(df[!, colname])))
    return df
end