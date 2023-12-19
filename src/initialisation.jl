"""
    initialize_population(::Type{T}, ::Type{S}, input; rng=Random.GLOBAL_RNG) where {T<: AbstractCell, S<: ModuleStructure}

Create the initial `Population` or `SinglelevelPopulation`. 
"""
function initialize_population(
    ::Type{T}, 
    ::Type{S}, 
    input::MultilevelInput;
    rng=Random.GLOBAL_RNG
) where {T<: AbstractCell, S<: ModuleStructure}

    N = getNinit(input)
    modulesize = getmaxmodulesize(input)
    Nmodules = getNmodules_init(input)
    modules = moduletype(T,S)[initialize(T, S, input.clonalmutations, N; rng) for _ in 1:Nmodules]
    homeostatic_modules, growing_modules = 
        if N == modulesize
            modules, moduletype(T,S)[]
        else
            moduletype(T,S)[], modules
        end

        return Population(
            homeostatic_modules, 
            growing_modules, 
            input.birthrate, input.deathrate, input.moranrate, input.asymmetricrate
        )
end

function initialize_population(
    ::Type{T}, 
    ::Type{S}, 
    input::SinglelevelInput;
    rng=Random.GLOBAL_RNG
) where {T<: AbstractCell, S<: ModuleStructure}

    N = getNinit(input)
    singlemodule = initialize(T, S, input.clonalmutations, N; rng)
    birthrate, deathrate, moranrate, asymmetricrate = getinputrates(input)

    return SinglelevelPopulation(singlemodule, birthrate, deathrate, moranrate, asymmetricrate)
end



"""
    initialize(::Type{T}, ::Type{S}, input, rng::AbstractRNG=Random.GLOBAL_RNG)

Initialise population of cells based on `input` and return as a `TreeModule{T, S} or
if `T == Cell` as a`CellModule{S}`.
"""
function initialize(
    ::Type{T}, 
    ::Type{S},
    clonalmutations,
    N;
    rng::AbstractRNG=Random.GLOBAL_RNG,
) where {T <: AbstractCell, S <: ModuleStructure}

    modulestructure = create_modulestructure(S, N)
    cells = create_cells(T, modulestructure, clonalmutations, N; rng)
    return new_module_from_cells(
        cells,
        0.0,
        Float64[0.0],
        1,
        0,
        modulestructure
    )
end

getNinit(input::Union{BranchingInput, BranchingMoranInput}) = 1
getNinit(input::MoranInput) = input.N
getNinit(input::MultilevelInput) = 1

getNmodules_init(input::MultilevelMoranInput) = input.maxmodules
getNmodules_init(input::Union{MultilevelBranchingInput, MultilevelBranchingMoranInput}) = 1

getmaxmodulesize(input::MultilevelInput) = input.modulesize
getmaxmodulesize(input::Union{MoranInput, BranchingMoranInput}) = input.N
getmaxmodulesize(input::BranchingInput) = Inf


function newcell(::Type{Cell}, id, mutations)
    return Cell(
        collect(1:mutations), 
        1,  #clonetype (wild-type)
        0,  #birthtime
        id, #unique cell id
        0   #parent id
    )
end

function newcell(::Type{T}, id, mutations) where T <: AbstractTreeCell
    return BinaryNode{T}(T(;id, mutations))
end

function create_cells(::Type{T}, structure::ModuleStructure, initialmutations, N=1; 
    rng=Random.GLOBAL_RNG) where T <:AbstractCell

    alivecells = [newcell(T, id, initialmutations) for id in 1:N]
    return position_cells(alivecells, structure, rng)
end

position_cells(alivecells, structure, rng) = alivecells

function position_cells(cells, structure::Linear, rng)
    N = length(cells)
    pad1 = round(Int64, (structure.size - N - 1)/2)
    pad2 = structure.size - pad1 - 1
    if rand(rng, 1:2) == 2 
        pad1, pad2 = pad2, pad1 
    end
    return [fill(nothing, pad1); cells; fill(nothing, pad2)]
end

function new_module_from_cells(cells::T, t, branchtimes, id, parentid, modulestructure::S) where {T, S}
    cellmodule = moduletype(T, S)(
        cells,
        t,
        branchtimes,
        id,
        parentid,
        modulestructure
    )
    return cellmodule
end