function getVAFresult(simulation::Simulation, rng::AbstractRNG=Random.GLOBAL_RNG; read_depth=100.0, 
    detectionlimit=5/read_depth, cellularity=1.0)

    trueVAF = getallelefreq(simulation)
    sampledVAF = sampledallelefreq(trueVAF, rng, read_depth=read_depth, 
        detectionlimit=detectionlimit, cellularity=cellularity)
    freq, freqp = subclonefreq(simulation.output)

    return VAFResult(
        read_depth,
        cellularity,
        detectionlimit,
        simulation.input,
        trueVAF,
        sampledVAF,
        freq,
        freqp
    )
end

function getVAFresult(multisim::MultiSimulation, rng::AbstractRNG=Random.GLOBAL_RNG; read_depth=100.0, 
    detectionlimit=5/read_depth, cellularity=1.0)

    trueVAFs = Array{Float64, 1}[]
    sampledVAFs = Array{Float64, 1}[]
    freqs = Array{Float64, 1}[]
    freqps = Array{Float64, 1}[]
    
    for moduletracker in multisim
        trueVAF = getallelefreq(moduletracker)
        sampledVAF = sampledallelefreq(trueVAF, rng, read_depth=read_depth, 
            detectionlimit=detectionlimit, cellularity=cellularity)
        freq, freqp = subclonefreq(moduletracker)
        push!(trueVAFs, trueVAF)
        push!(sampledVAFs, sampledVAF)
        push!(freqs, freq)
        push!(freqps, freqp)
    end

    return VAFResultMulti(
        read_depth,
        cellularity,
        detectionlimit,
        simulation.input,
        trueVAFs,
        sampledVAFs,
        freqs,
        freqps
    )
end

function getallelefreq(simulation::Simulation)
    return getallelefreq(simulation.output, simulation.input.ploidy)
end

function getallelefreq(moduletracker::ModuleTracker, ploidy)
    mutations = cellsconvert(moduletracker.cells).mutations
    return getallelefreq(mutations, moduletracker.Nvec[end], ploidy)
end

function getallelefreq(mutations, N, ploidy=2)
    allelefreq = counts(mutations,minimum(mutations):maximum(mutations))
    # idx = f .> 0.01 #should this be f/(2*cellnum) .> 0.01, i.e. only include freq > 1% ??
    allelefreq = map(Float64, allelefreq)
    allelefreq ./= (ploidy * N) #correct for ploidy
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

"""
    sampledhist(trueVAF::Array{Float64,1}, cellnum::Int64; <keyword arguments>)

Create synthetic experimental data by smapling from the true VAF distribution 
according to experimental constraints. 

### Keyword arguments
- `detectionlimit = 5/read_depth`: minimum VAF which can be detected
- `read_depth = 100`: expected number of reads for each allele
- `cellularity = 1.0`: 

"""
function sampledallelefreq(trueVAF::Array{Float64, 1}, rng::AbstractRNG; 
    read_depth=100.0, detectionlimit=5/read_depth, cellularity=1.0)

    VAF = trueVAF * cellularity
    filter!(x -> x > detectionlimit, VAF)
    # Pois(np) ≈ B(n,p) when n >> p, thus we can sample from Poisson(read_depth)
    # rather than Binomial(cellnum, read_depth/cellnum), length(VAF)) as D << N^2
    depth = rand(rng, Poisson(read_depth), length(VAF))
    sampalleles = map((n, p) -> rand(rng, Binomial(n, p)), depth, VAF)
    sampledVAF = sampalleles ./ depth
    return sampledVAF
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

function subclonefreq(moduletracker)
    #get proportion of cells in each subclone
    clonesize = getclonesize(moduletracker)
    clonefreqp = clonesize[2:end]/sum(clonesize)
    clonefreq = copy(clonefreqp)
    if length(clonefreqp) > 1
        clonefreq, subclonalmutations = 
            calculateclonefreq!(clonefreq, subclonalmutations, moduletracker.subclones)
    end
    return clonefreq, clonefreqp
end