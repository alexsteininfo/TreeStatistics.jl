function getVAFresult(simulation, rng::AbstractRNG=Random.GLOBAL_RNG; read_depth=100.0, 
    detectionlimit=5/read_depth, cellularity=1.0)
    trueVAF = getallelefreq(simulation)
    sampledVAF = sampledallelefreq(trueVAF, rng, read_depth=read_depth, 
        detectionlimit=detectionlimit, cellularity=cellularity)
    freq, freqp = subclonefreq(simulation.output.subclones)

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

function getVAFresult(multisimulation::MultiSimulation, moduleid, rng::AbstractRNG=Random.GLOBAL_RNG; read_depth=100.0, 
    detectionlimit=5/read_depth, cellularity=1.0)

    trueVAF = getallelefreq(multisimulation, moduleid)
    sampledVAF = sampledallelefreq(trueVAF, rng, read_depth=read_depth, 
        detectionlimit=detectionlimit, cellularity=cellularity)
    freq, freqp = subclonefreq(multisimulation[moduleid], length(multisimulation.output.subclones))

    return VAFResult(
        read_depth,
        cellularity,
        detectionlimit,
        multisimulation.input,
        trueVAF,
        sampledVAF,
        freq,
        freqp
    )
end

function getVAFresult(multisimulation::MultiSimulation, rng::AbstractRNG=Random.GLOBAL_RNG; read_depth=100.0, 
    detectionlimit=5/read_depth, cellularity=1.0)

    trueVAF = getallelefreq(multisimulation)
    sampledVAF = sampledallelefreq(trueVAF, rng, read_depth=read_depth, 
        detectionlimit=detectionlimit, cellularity=cellularity)
    freq, freqp = Float64[0.0], Float64[0.0] #TODO not implemented properly

    return VAFResult(
        read_depth,
        cellularity,
        detectionlimit,
        multisimulation.input,
        trueVAF,
        sampledVAF,
        freq,
        freqp
    )
end

function getVAFresultmulti(multisim::MultiSimulation, rng::AbstractRNG=Random.GLOBAL_RNG; read_depth=100.0, 
    detectionlimit=5/read_depth, cellularity=1.0, ploidy=2)

    trueVAFs = Vector{Float64}[]
    sampledVAFs = Vector{Float64}[]
    freqs = Vector{Float64}[]
    freqps = Vector{Float64}[]
    
    for cellmodule in multisim
        trueVAF = getallelefreq(cellmodule, multisim.input.ploidy)
        sampledVAF = sampledallelefreq(trueVAF, rng, read_depth=read_depth, 
            detectionlimit=detectionlimit, cellularity=cellularity)
        freq, freqp = subclonefreq(cellmodule)
        push!(trueVAFs, trueVAF)
        push!(sampledVAFs, sampledVAF)
        push!(freqs, freq)
        push!(freqps, freqp)
    end

    return VAFResultMulti(
        read_depth,
        cellularity,
        detectionlimit,
        multisim.input,
        trueVAFs,
        sampledVAFs,
        freqs,
        freqps
    )
end

"""
    getallelefreq(simulation)

Calculate allele frequency for each mutation. If `simulation <: MultiSimulation` the allele
    freq is calculated for all alleles in all modules.
"""
function getallelefreq(simulation)
    return getallelefreq(simulation.output, simulation.input.ploidy)
end

"""
    getallelefreq(population::MultiSimulation, moduleid)

Calculate allele frequency for each mutation in the module (or vector of modules) specified 
    by `moduleid`.
"""
function getallelefreq(simulation::Simulation, moduleid)
    return getallelefreq(simulation.output[moduleid], simulation.input.ploidy)
end

"""
    getallelefreq(module::AbstractModule, ploidy)
    getallelefreq(cellvector::AbstractCellVector, ploidy)


Calculate allele frequency for all mutations in a single `module` or `cellvector` with given
    `ploidy`.
"""
function getallelefreq(singlelevelpopulation::SinglelevelPopulation, ploidy)
    return getallelefreq(singlelevelpopulation.singlemodule.cells, ploidy)
end

function getallelefreq(cells::CellVector, ploidy)
    N = length(cells)
    mutations = cellsconvert(cells).mutations
    allelefreq = length(mutations) == 0 ? Float64[] : Float64.(counts(mutations))
    filter!(x -> x != 0, allelefreq)
    allelefreq ./= (ploidy * N) #correct for ploidy
    return allelefreq
end

function getallelefreq(cells::AbstractTreeCellVector, ploidy)
    N = length(cells)
    allelefreqs = Float64[]
    nodes = getroot(cells)
    node_freqs = Float64[cell_subset_size(node, cells) / (ploidy*N) for node in nodes]
    while true
        #get next node and freq
        node = popfirst!(nodes) 
        node_freq = popfirst!(node_freqs)
        for i in 1:node.data.mutations
            push!(allelefreqs, node_freq)
        end
        #add any live child nodes to nodes list
        left_alive, right_alive = isalive(node.left), isalive(node.right)
        if left_alive
            if right_alive
                n_descendents_left = cell_subset_size(node.left, cells)
                if n_descendents_left != 0
                    push!(nodes, node.left)
                    push!(node_freqs, n_descendents_left / (ploidy*N))
                end
                n_descendents_right = cell_subset_size(node.right, cells)
                if n_descendents_right != 0
                    push!(nodes, node.right)
                    push!(node_freqs, n_descendents_right / (ploidy*N))
                end
            else
                push!(nodes, node.left)
                push!(node_freqs, node_freq) #if only one child freq is same as for parent
            end
        else
            if right_alive
                push!(nodes, node.right)
                push!(node_freqs, node_freq) #if only one child freq is same as for parent
            else
                #if no new nodes have been added check whether any nodes are left
                length(nodes) == 0 && break 
            end
        end
    end
    return allelefreqs
end

"""
    getallelefreq(module::Vector{AbstractModule}, ploidy)

Calculate allele frequency for all mutations in a vector of `modules` with given `ploidy`.
"""
function getallelefreq(population::Population, ploidy)
    return getallelefreq([cell for m in population for cell in m.cells], ploidy)
end

function getfixedallelefreq(mutations::Vector{Int64})
    return counts(mutations)    
end

getfixedallelefreq(simulation, idx=nothing) = getfixedallelefreq(simulation.output, idx)

function getfixedallelefreq(population::Population{T}, idx=nothing) where T <: CellModule
    mutations = reduce(vcat, clonal_mutation_ids(population, idx))
    return getfixedallelefreq(mutations)
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
    sampledallelefreq(trueVAF::Vector{Float64}, rng::AbstractRNG; 
    read_depth=100.0, detectionlimit=5/read_depth, cellularity=1.0)

Create synthetic experimental data by smapling from the true VAF distribution 
according to experimental constraints. 

### Keyword arguments
- `detectionlimit = 5/read_depth`: minimum VAF which can be detected
- `read_depth = 100`: expected number of reads for each allele
- `cellularity = 1.0`: 

"""
function sampledallelefreq(trueVAF::Vector{Float64}, rng::AbstractRNG; 
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
    gethist(VAF::Vector{Float64}; fmin = 0.0, fmax = 1.0, fstep = 0.001)

Fit VAF data into bins `fmin`:`fstep`:`fmax` and return DataFrame with columns
:VAF, :freq [= ``m(f)``], :cumfreq [= ``M(f)``].

Note ``M(f)`` is the cummulative number of mutations with frequency ``f``, i.e. the 
number of mutations with frequency in ``(f, 1)``, while ``m(f)`` is the number of mutations with 
frequency in ``(f - fstep, f)``.

"""
function gethist(VAF::Vector{Float64}; fmin = 0.0, fmax = 1, fstep = 0.001, norm=:density)
    x = fmin:fstep:fmax
    y = fit(Histogram, VAF, x; closed=:right).weights
    if norm == :density
        y ./= fstep
    end
    dfhist = DataFrame(VAF = x[2:end], freq = y)
    dfhist = addcumfreq!(dfhist,:freq)
    return dfhist
end


function addcumfreq!(df,colname)
    df[!, Symbol(:cum,colname)] = reverse(cumsum(reverse(df[!, colname])))
    return df
end


function subclonefreq(cellmodule, nsubclones)
    #get proportion of cells in each subclone
    clonesize = getsubclonesizes(cellmodule, nsubclones)
    clonefreqp = clonesize[2:end]/sum(clonesize)
    clonefreq = copy(clonefreqp)
    if length(clonefreqp) > 1
        clonefreq, subclonalmutations = 
            calculateclonefreq!(clonefreq, subclonalmutations, cellmodule.subclones)
    end
    return clonefreq, clonefreqp
end

function subclonefreq(population::Vector)
    #TODO need to implement this
    Float64[], Float64[]
end

function subclonefreq(treemodule::TreeModule)
    #selection not implemented for tree type
    return Float64[], Float64[]
end