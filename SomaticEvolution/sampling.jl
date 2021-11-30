
"""
    sampledhist(trueVAF::Array{Float64,1}, cellnum::Int64; <keyword arguments>)

Create synthetic experimental data by smapling from the true VAF distribution 
according to experimental constraints. 

### Keyword arguments
- `detectionlimit = 5/read_depth`: minimum VAF which can be detected
- `read_depth = 100`: expected number of reads for each allele
- `cellularity = 1.0`: 

"""
function sampledhist(trueVAF::Array{Float64, 1}, cellnum::Int64, rng::AbstractRNG; 
    detectionlimit=5/read_depth, read_depth=100.0, cellularity=1.0)

    VAF = trueVAF * cellularity
    filter!(x -> x > detectionlimit, VAF)
    # Pois(np) ≈ B(n,p) when n >> p, thus we can sample from Poisson(read_depth)
    # rather than Binomial(cellnum, read_depth/cellnum), length(VAF)) as D << N^2
    depth = rand(rng, Poisson(read_depth), length(VAF))
    sampalleles = map((n, p) -> rand(rng, Binomial(n, p)), depth, VAF)
    sampledVAF = sampalleles ./ depth
    return SampledData(sampledVAF, depth)
end

"""
    gethist(VAF::Array{Float64,1}; fmin = 0.0, fmax = 1.0, fstep = 0.001)

Fit VAF data into bins `fmin`:`fstep`:`fmax` and return DataFrame with columns
:VAF, :freq [= ``m(f)``], :cumfreq [= ``M(f)``].

Note ``M(f)`` is the cummulative number of mutations with frequency ``f``, i.e. the 
number of mutations with frequency in ``(f, 1)``, while ``m(f)`` is the number of mutations with 
frequency in ``(f - fstep, f)``.

"""
function gethist(VAF::Array{Float64,1}; fmin = 0.0, fmax = 1, fstep = 0.001)
    x = fmin:fstep:fmax
    y = fit(Histogram, VAF, x, closed=:right).weights
    dfhist = DataFrame(VAF = x[2:end], freq = y)
    dfhist = addcumfreq!(dfhist,:freq)
    return dfhist
end


function addcumfreq!(df,colname)
    df[!, Symbol(:cum,colname)] = reverse(cumsum(reverse(df[!, colname])))
    return df
end