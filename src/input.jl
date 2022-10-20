"""
    SimulationInput

Supertype for simulation inputs.

# Subtypes:
    - `BranchingInput <: SinglelevelInput`: input for simulating a branching process
    - `MoranInput <: SinglelevelInput`: input for simulating a Moran process 
    - `BranchingMoranInput <: SinglelevelInput`: input for simulating a branching process 
        until fixed size is reached, then switching to Moran process
    - `MultilevelBranchingInput <: MultilevelInput`: input for simulating a module branching 
        process
    - `MultilevelBranchingMoranInput <: MultilevelInput` input for simulating a module 
        branching process until fixed size is reached, then switching to Moran
"""
abstract type SimulationInput end

abstract type SinglelevelInput <: SimulationInput end

"""
    BranchingInput <: SinglelevelInput <: SimulationInput
"""
struct BranchingInput <: SinglelevelInput
    numclones::Int64 
    Nmax::Int64
    tmax::Float64
    clonalmutations::Int64
    selection::Vector{Float64}
    μ::Float64
    b::Float64
    d::Float64
    tevent::Vector{Float64}
    mutationdist::Symbol
    maxclonesize::Union{Int64, Nothing}
    ploidy::Int64
end

"""
    MoranInput <: SinglelevelInput <: SimulationInput
"""
struct MoranInput <: SinglelevelInput
    N::Int64
    numclones::Int64 
    tmax::Float64
    clonalmutations::Int64
    selection::Vector{Float64}
    μ::Float64
    bdrate::Float64
    tevent::Vector{Float64}
    mutationdist::Symbol
    ploidy::Int64
end

"""
    BranchingMoranInput <: SinglelevelInput <: SimulationInput
"""
struct BranchingMoranInput <: SinglelevelInput
    numclones::Int64 
    Nmax::Int64
    tmax::Float64
    clonalmutations::Int64
    selection::Vector{Float64}
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    tevent::Vector{Float64}
    mutationdist::Symbol
    ploidy::Int64
end

abstract type MultilevelInput <: SimulationInput end

"""
    MultilevelBranchingInput <: MultilevelInput <: SimulationInput
"""
struct MultilevelBranchingInput <: MultilevelInput
    modulesize::Int64
    tmax::Float64
    maxmodules::Int64
    clonalmutations::Int64
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    mutationdist::Symbol
    branchrate::Float64
    branchinitsize::Int64
    ploidy::Int64
end

"""
    MultilevelBranchingMoranInput <: MultilevelInput <: SimulationInput
"""
struct MultilevelBranchingMoranInput <: MultilevelInput
    modulesize::Int64
    tmax::Float64
    maxmodules::Int64
    clonalmutations::Int64
    μ::Float64
    bdrate::Float64
    b::Float64
    d::Float64
    mutationdist::Symbol
    branchrate::Float64
    branchinitsize::Int64
    ploidy::Int64
end

"""
    BranchingInput(<keyword arguments>)

Input for a single level branching process simulation that starts with a single cell.
    
# Keyword arguments:
- `Nmax::Int64 = 10000`: maximum number of cells
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Float64 = 10.0`: mutation rate per division per cell
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `b::Float64 = log(2.0)`: birth rate for wild-type cells
- `d::Float64 = 0.0`: death rate for wild-type cells
- `numclones::Int64 = 0`: number of mutant subclones
- `selection::Vector{Float64} = [0.0, 0.0, ...]`: selection strength of each mutant subclone
- `tevent::Vector{Float64} = [1.0, 1.5, ...]`: time each mutant arises
- `fixedmu::Bool = false`: if `mutantdist` is not specified defines whether the number of 
    mutations accumulated at division is :poisson or :fixed
- `mutationdist::Symbol = fixedmu ? :fixed : :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `maxclonesize::Int64 = nothing`: ...
"""
function BranchingInput(;numclones=0, Nmax=10000, tmax=Inf, ploidy=2, μ=10.0, 
    clonalmutations=0, selection=fill(0.0,numclones), b=log(2.0), d=0.0, 
    tevent=collect(1.0:0.5:(1+numclones)/2), fixedmu=false, mutationdist=nothing,
    maxclonesize=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

    return BranchingInput(
        numclones,
        Nmax,
        tmax,
        clonalmutations,
        selection,
        μ,
        b,
        d,
        tevent,
        mutationdist,
        maxclonesize,
        ploidy
    )
end

"""
    MoranInput(<keyword arguments>)

Input for a single level Moran process simulation that starts with a `N` identical cells.

# Keyword arguments:
- `N::Int64 = 10000`: number of cells
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Float64 = 10.0`: mutation rate per division per cell
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `bdrate::Float64 = log(2.0)`: birth/death rate for wild-type cells
- `numclones::Int64 = 0`: number of mutant subclones
- `selection::Vector{Float64} = [0.0, 0.0, ...]`: selection strength of each mutant subclone
- `tevent::Vector{Float64} = [1.0, 1.5, ...]`: time each mutant arises
- `fixedmu::Bool = false`: if `mutantdist` is not specified defines whether the number of 
    mutations accumulated at division is :poisson or :fixed
- `mutationdist::Symbol = fixedmu ? :fixed : :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
"""
function MoranInput(;numclones=0, N=10000, ploidy=2, μ=10.0, clonalmutations=0, 
    selection=fill(0.0,numclones), bdrate=log(2.0), tmax=15.0,
    tevent=collect(1.0:0.5:(1+numclones)/2), fixedmu=false, mutationdist=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)

    return MoranInput(
        N,
        numclones,
        tmax,
        clonalmutations,
        selection,
        μ,
        bdrate,
        tevent,
        mutationdist,
        ploidy
    )
end

"""
    BranchingMoranInput(<keyword arguments>)

Input for a single level simulation that starts with a single cell and simulated a 
    branching process until the population reaches `Nmax`, then switches to a Moran process.

# Keyword arguments:
- `Nmax::Int64 = 10000`: maximum number of cells
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Float64 = 10.0`: mutation rate per division per cell
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `b::Float64 = log(2.0)`: birth rate for wild-type cells in branching phase
- `d::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `bdrate::Float64 = log(2.0)`: birth/death rate for wild-type cells in Moran phase
- `numclones::Int64 = 0`: number of mutant subclones
- `selection::Vector{Float64} = [0.0, 0.0, ...]`: selection strength of each mutant subclone
- `tevent::Vector{Float64} = [1.0, 1.5, ...]`: time each mutant arises
- `fixedmu::Bool = false`: if `mutantdist` is not specified defines whether the number of 
    mutations accumulated at division is :poisson or :fixed
- `mutationdist::Symbol = fixedmu ? :fixed : :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
"""
function BranchingMoranInput(;numclones=1, Nmax=10000, ploidy=2, μ=10.0, 
    clonalmutations=0, selection=fill(0.0,numclones), bdrate=nothing, b=nothing, 
    d=0, tmax=15.0, tevent=collect(1.0:0.5:(1+numclones)/2), fixedmu=false, 
    mutationdist=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)
    b, bdrate = set_cell_birthrates(b, bdrate)

    return BranchingMoranInput(
        numclones,
        Nmax,
        tmax,
        clonalmutations,
        selection,
        μ,
        bdrate,
        b,
        d,
        tevent,
        mutationdist,
        ploidy
    )
    
end

"""
    MultilevelBranchingInput(<keyword arguments>)

Input for a multilevel branching simulation that starts with a single cell in a single module. 
    
Within module dynamics follows a branching process until `modulesize` is reached and then
switches to a Moran process. Module level dynamics follow a branching process (homeostatic 
modules branch at rate `branchrate`) with no death. 

# Keyword arguments:
- `modulesize::Int64 = 200`: maximum number of cells per module
- `maxmodules::Int64 = 10000`: maximum number of modules in population 
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Float64 = 10.0`: mutation rate per division per cell
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `b::Float64 = 1.0`: birth rate for wild-type cells in branching phase
- `d::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `bdrate::Float64 = 1.0`: birth/death rate for wild-type cells in Moran phase
- `branchrate::Float64 = 5.0`: rate at which homeostatic modules split to form new modules
- `branchfraction::Float64 = 0.1`: fraction of cells sampled to form a new module (only if
    `branchinitsize` is not given)
- `branchinitsize::Int64 = branchfraction * modulesize`: number of cells sampled to form a 
    new module
- `fixedmu::Bool = false`: if `mutantdist` is not specified defines whether the number of 
    mutations accumulated at division is :poisson or :fixed
- `mutationdist::Symbol = fixedmu ? :fixed : :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
"""
function MultilevelBranchingInput(;modulesize=200, ploidy=2, μ=10.0, clonalmutations=0, 
    bdrate=nothing, b=nothing, d=0, tmax=15, maxmodules=10000, fixedmu=false, 
    mutationdist=nothing, branchrate=5, branchfraction=0.1, branchinitsize=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)
    b, bdrate = set_cell_birthrates(b, bdrate)

    return MultilevelBranchingInput(
            modulesize,
            tmax,
            maxmodules,
            clonalmutations,
            μ,
            bdrate,
            b,
            d,
            mutationdist,
            branchrate,
            branchinitsize !== nothing ? branchinitsize : ceil(modulesize * branchfraction),
            ploidy
    )
end

"""
    MultilevelBranchingMoranInput(<keyword arguments>)

Input for a multilevel simulation of a homeostatic population that starts with a single cell 
    in a single module. 
    
Within module dynamics follows a branching process until `modulesize` is reached and then
switches to a Moran process. Module level dynamics follow a branching process (homeostatic 
modules branch at rate `branchrate`) with no death. Once module population reaches 
`maxmodules`, switch to a Moran process.

# Keyword arguments:
- `modulesize::Int64 = 200`: maximum number of cells per module
- `maxmodules::Int64 = 10000`: maximum number of modules in population 
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Float64 = 10.0`: mutation rate per division per cell
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `b::Float64 = 1.0`: birth rate for wild-type cells in branching phase
- `d::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `bdrate::Float64 = 1.0`: birth/death rate for wild-type cells in Moran phase
- `branchrate::Float64 = 5.0`: rate at which homeostatic modules split to form new modules
- `branchfraction::Float64 = 0.1`: fraction of cells sampled to form a new module (only if
    `branchinitsize` is not given)
- `branchinitsize::Int64 = branchfraction * modulesize`: number of cells sampled to form a 
    new module
- `fixedmu::Bool = false`: if `mutantdist` is not specified defines whether the number of 
    mutations accumulated at division is :poisson or :fixed
- `mutationdist::Symbol = fixedmu ? :fixed : :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
"""
function MultilevelBranchingMoranInput(;modulesize=200, ploidy=2, μ=10.0, clonalmutations=0, 
    bdrate=nothing, b=nothing, d=0, tmax=15, maxmodules=10000, fixedmu=false, 
    mutationdist=nothing, branchrate=0.1, branchfraction=0.1, branchinitsize=nothing)

    mutationdist = set_mutationdist(mutationdist, fixedmu)
    b, bdrate = set_cell_birthrates(b, bdrate)

    return MultilevelBranchingMoranInput(
            modulesize,
            tmax,
            maxmodules,
            clonalmutations,
            μ,
            bdrate,
            b,
            d,
            mutationdist,
            branchrate,
            branchinitsize !== nothing ? branchinitsize : ceil(modulesize * branchfraction),
            ploidy,
    )
end

"""
    newinput(::Type{InputType}; kwargs) where InputType <: SimulationInput

Create a new input of type `InputType`.
"""
function newinput(::Type{InputType}; kwargs) where InputType <: SimulationInput
    return InputType(;kwargs...)
end

"""
    newinput(input::InputType; kwargs...) where InputType <: SimulationInput

Create a new input of type `InputType`. Any fields not given in `kwargs` default to the 
values in `input`.
    """
function newinput(input::InputType; kwargs...) where InputType <: SimulationInput
    newkwargs = Dict(
        field in keys(kwargs) ? field => kwargs[field] : field => getfield(input, field)
            for field in fieldnames(InputType))
    return InputType(;newkwargs...)
end

function set_mutationdist(mutationdist, fixedmu)
    if isnothing(mutationdist)
        return fixedmu ? :fixed : :poisson
    elseif typeof(mutationdist) === String
        return Symbol(mutationdist)
    else 
        return mutationdist
    end
end

function set_cell_birthrates(b, bdrate)
    return begin
        if isnothing(b) && isnothing(bdrate)
            1.0, 1.0
        elseif isnothing(b)
            bdrate, bdrate
        elseif isnothing(bdrate)
            b, b
        else
            b, bdrate
        end
    end
end