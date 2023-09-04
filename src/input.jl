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

#region Single-level simulation inputs
abstract type SinglelevelInput <: SimulationInput end

"""
    BranchingInput <: SinglelevelInput <:SimulationInput

Input for a single level branching process simulation that starts with a single cell.
    
# Keyword arguments:
- `Nmax::Int64 = 1000`: maximum number of cells
- `tmax::Float64 = Inf`: maximum time to run simulation
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells
- `deathrate::Float64 = 0.0`: death rate for wild-type cells
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Float64 = 1.0`: mutation rate per division per cell
- `mutationdist::Symbol = :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `ploidy::Int64 = 2`: cell ploidy 
- `numclones::Int64 = 0`: number of mutant subclones
- `selection::Vector{Float64} = [0.0, 0.0, ...]`: selection strength of each mutant subclone
- `tevent::Vector{Float64} = [1.0, 1.5, ...]`: time each mutant arises
"""
Base.@kwdef struct BranchingInput <: SinglelevelInput
    Nmax::Int64 = 1000
    tmax::Float64 = Inf
    birthrate::Float64 = 1.0
    deathrate::Float64 = 0.0
    clonalmutations::Int64 = 0
    μ::Float64 = 1.0
    mutationdist::Symbol = :poisson
    ploidy::Int64 = 2
    numclones::Int64 = 0
    selection::Vector{Float64} = fill(0.0,numclones)
    tevent::Vector{Float64} = collect(1.0:0.5:(1+numclones)/2)
end

"""
    MoranInput <: SinglelevelInput

Input for a single level Moran process simulation that starts with `N` identical cells.

# Keyword arguments:
- `N::Int64 = 1000`: number of cells
- `tmax::Float64 = 15.0`: maximum time to run simulation
- `moranrate::Float64 = 1.0`: moran update rate for wild-type cells
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Float64 = 1.0`: mutation rate per division per cell
- `mutationdist::Symbol = :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `ploidy::Int64 = 2`: cell ploidy 
- `numclones::Int64 = 0`: number of mutant subclones
- `selection::Vector{Float64} = [0.0, 0.0, ...]`: selection strength of each mutant subclone
- `tevent::Vector{Float64} = [1.0, 1.5, ...]`: time each mutant arises
"""
Base.@kwdef struct MoranInput <: SinglelevelInput
    N::Int64 = 1000
    tmax::Float64 = 15.0
    moranrate::Float64 = 1.0
    moranincludeself::Bool = true
    clonalmutations::Int64 = 0
    μ::Float64 = 1.0
    mutationdist::Symbol = :poisson
    ploidy::Int64 = 2
    numclones::Int64 = 0
    selection::Vector{Float64} = fill(0.0,numclones)
    tevent::Vector{Float64} = collect(1.0:0.5:(1+numclones)/2)
end

"""
    BranchingMoranInput <: SinglelevelInput

Input for a single level simulation that grows by a branching process to `Nmax` cells and
    then switches to a Moran process.

# Keyword arguments:
- `N::Int64 = 1000`: number of cells
- `tmax::Float64 = 15.0`: maximum time to run simulation
- `moranrate::Float64 = 1.0`: moran update rate for wild-type cells
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `birthrate::Float64 = moranrate`: birth rate for wild-type cells
- `deathrate::Float64 = 0.0`: death rate for wild-type cells
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Float64 = 1.0`: mutation rate per division per cell
- `mutationdist::Symbol = :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `ploidy::Int64 = 2`: cell ploidy
- `numclones::Int64 = 0`: number of mutant subclones
- `selection::Vector{Float64} = [0.0, 0.0, ...]`: selection strength of each mutant subclone
- `tevent::Vector{Float64} = [1.0, 1.5, ...]`: time each mutant arises
"""
Base.@kwdef struct BranchingMoranInput <: SinglelevelInput
    Nmax::Int64 = 1000
    tmax::Float64 = 15.0
    moranrate::Float64 = 1.0
    moranincludeself::Bool = true
    birthrate::Float64 = moranrate
    deathrate::Float64 = 0.0
    clonalmutations::Int64 = 0
    μ::Float64 = 1.0
    mutationdist::Symbol = :poisson
    ploidy::Int64 = 2
    numclones::Int64 = 0
    selection::Vector{Float64} = fill(0.0,numclones)
    tevent::Vector{Float64} = collect(1.0:0.5:(1+numclones)/2)
end
#endregion

#region Multi-level simulation inputs
abstract type MultilevelInput <: SimulationInput end

"""
    MultilevelBranchingInput <: MultilevelInput

Input for a multilevel branching simulation that starts with a single cell in a single module. 
    
Within module dynamics follows a branching process until `modulesize` is reached and then
switches to a Moran process. Module level dynamics follow a branching process (homeostatic 
modules branch at rate `branchrate`) with no death. 

# Keyword arguments:
- `modulesize::Int64 = 20`: maximum number of cells per module
- `maxmodules::Int64 = 1000`: maximum number of modules in population 
- `tmax::Float64 = Inf`: maximum time to run simulation
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells in branching phase
- `deathrate::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `moranrate::Float64 = 1.0`: rate of Moran updating for wild-type cells in homeostasis
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `asymmetricrate::Float64 = 0.0`: rate of asymmetric updating for wild-type cells in homeostasis
- `branchrate::Float64 = 5.0`: rate at which homeostatic modules split to form new modules
- `branchinitsize::Int64 = 1`: number of cells sampled to form a 
    new module
- `modulebranching::Symbol = :split`: determines the method by which a new module is formed
    at branching. Options are `:split` (module cells are split between two modules), 
    `:withreplacement` (cells are sampled from the parent module and undergo division with one
    cell returning to the parent, before the next cell is sampled, and the other entering 
    the new module), `:withoutreplacement` (as previous except that cells are returned to
    parent module after all smapling is completed), `:withreplacement_nomutations` and
    `withoutreplacement_nomutations` (as previous but dividing cells get no new mutations).
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Float64 = 1.0`: mutation rate per division per cell
- `mutationdist::Symbol = :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `ploidy::Int64 = 2`
"""
Base.@kwdef struct MultilevelBranchingInput <: MultilevelInput
    modulesize::Int64 = 20
    maxmodules::Int64 = 1000
    tmax::Float64 = Inf
    moranrate::Float64 = 1.0
    moranincludeself::Bool = true
    asymmetricrate::Float64 = 0.0
    birthrate::Float64 = 1.0
    deathrate::Float64 = 0.0
    branchrate::Float64  = 5.0
    branchinitsize::Int64 = 1
    modulebranching::Symbol = :split
    clonalmutations::Int64 = 0
    μ::Float64 = 1.0
    mutationdist::Symbol = :poisson
    ploidy::Int64 = 2
end

"""
    MultilevelMoranInput <: MultilevelInput

Input for a multilevel branching simulation that starts with `maxmodules` modules, each with
a single cell. 
    
Within module dynamics follows a branching process until `modulesize` is reached and then
switches to a Moran process. Module level dynamics follows a Moran process at rate 
`branchrate`.

# Keyword arguments:
- `modulesize::Int64 = 200`: maximum number of cells per module
- `maxmodules::Int64 = 10000`: maximum number of modules in population 
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Float64 = 10.0`: mutation rate per division per cell
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells in branching phase
- `deathrate::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `moranrate::Float64 = 1.0`: rate of Moran updating for wild-type cells in homeostasis
- `asymmetricrate::Float64 = 0.0`: rate of asymmetric updating for wild-type cells in homeostasis
- `branchrate::Float64 = 5.0`: rate at which homeostatic modules split to form new modules
- `branchinitsize::Int64 = 1`: number of cells sampled to form a 
    new module
- `mutationdist::Symbol = :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `modulebranching::Symbol = :split`: determines the method by which a new module is formed
    at branching. Options are `:split` (module cells are split between two modules), 
    `:withreplacement` (cells are sampled from the parent module and undergo division with one
    cell returning to the parent, before the next cell is sampled, and the other entering 
    the new module), `:withoutreplacement` (as previous except that cells are returned to
    parent module after all smapling is completed), `:withreplacement_nomutations` and
    `withoutreplacement_nomutations` (as previous but dividing cells get no new mutations).
"""


Base.@kwdef struct MultilevelMoranInput <: MultilevelInput
    modulesize::Int64 = 20
    maxmodules::Int64 = 1000
    tmax::Float64 = 15.0
    moranrate::Float64 = 1.0
    moranincludeself::Bool = true
    asymmetricrate::Float64 = 0.0
    birthrate::Float64 = 1.0
    deathrate::Float64 = 0.0
    branchrate::Float64  = 5.0
    branchinitsize::Int64 = 1
    modulebranching::Symbol = :split
    clonalmutations::Int64 = 0
    μ::Float64 = 1.0
    mutationdist::Symbol = :poisson
    ploidy::Int64 = 2
end


"""
    MultilevelBranchingMoranInput <: MultilevelInput

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
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells in branching phase
- `deathrate::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `moranrate::Float64 = 1.0`: birth/death rate for wild-type cells in Moran phase
- `asymmetricrate::Float64 = 0.0`: rate of asymmetric updating for wild-type cells in homeostasis
- `branchrate::Float64 = 5.0`: rate at which homeostatic modules split to form new modules
- `branchinitsize::Int64 = 1`: number of cells sampled to form a 
    new module
- `mutationdist::Symbol = :poisson`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric)
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `modulebranching::Symbol = :split`: determines the method by which a new module is formed
    at branching. Options are `:split` (module cells are split between two modules), 
    `:withreplacement` (cells are sampled from the parent module and undergo division with one
    cell returning to the parent, before the next cell is sampled, and the other entering 
    the new module), `:withoutreplacement` (as previous except that cells are returned to
    parent module after all smapling is completed), `:withreplacement_nomutations` and
    `withoutreplacement_nomutations` (as previous but dividing cells get no new mutations).
"""
Base.@kwdef struct MultilevelBranchingMoranInput <: MultilevelInput
    modulesize::Int64 = 20
    maxmodules::Int64 = 1000
    tmax::Float64 = 15.0
    moranrate::Float64 = 1.0
    moranincludeself::Bool = true
    asymmetricrate::Float64 = 0.0
    birthrate::Float64 = 1.0
    deathrate::Float64 = 0.0
    branchrate::Float64  = 5.0
    branchinitsize::Int64 = 1
    modulebranching::Symbol = :split
    clonalmutations::Int64 = 0
    μ::Float64 = 1.0
    mutationdist::Symbol = :poisson
    ploidy::Int64 = 2
end

#endregion


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
