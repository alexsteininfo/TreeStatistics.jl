"""
    SimulationInput

Supertype for simulation inputs.

# Subtypes
- `BranchingInput <: SinglelevelInput`: input for simulating a branching process
- `MoranInput <: SinglelevelInput`: input for simulating a Moran process 
- `BranchingMoranInput <: SinglelevelInput`: input for simulating a branching process 
    until fixed size is reached, then switching to Moran process
- `MultilevelBranchingInput <: MultilevelInput`: input for simulating a module branching 
    process (each module is formed of cells and grows by a branching process then
    switches to a Moran process and/or asymmetric cell division)
- `MultilevelBranchingMoranInput <: MultilevelInput` input for simulating a module 
    branching process until fixed size is reached, then switching to Moran (module-level
    dynamics are the same as for MultilevelBranchingInput)
"""
abstract type SimulationInput end

#region Single-level simulation inputs
abstract type SinglelevelInput <: SimulationInput end

"""
    BranchingInput <: SinglelevelInput <:SimulationInput

Input type for a single level branching process simulation that starts with a single cell.
    
# Keyword arguments:
- `Nmax::Int64 = 1000`: maximum number of cells
- `tmax::Float64 = Inf`: maximum time to run simulation
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells
- `deathrate::Float64 = 0.0`: death rate for wild-type cells
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes. 
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should 
    match `length(μ)`.
- `ploidy::Int64 = 2`: cell ploidy 
"""
struct BranchingInput <: SinglelevelInput
    Nmax::Int64
    tmax::Float64
    birthrate::Float64
    deathrate::Float64
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function BranchingInput(;
    Nmax = 1000,
    tmax = Inf,
    birthrate = 1.0,
    deathrate = 0.0,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)   
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"
    return BranchingInput(
        Nmax, tmax, birthrate, deathrate, clonalmutations, μ, mutationdist, ploidy
    )
end

tovector(a::T) where T = T[a]
tovector(a::Vector{T}) where T = a

"""
    MoranInput <: SinglelevelInput

Input type for a single level Moran process simulation that starts with `N` identical cells.

# Keyword arguments:
- `N::Int64 = 1000`: number of cells
- `tmax::Float64 = 15.0`: maximum time to run simulation
- `moranrate::Float64 = 1.0`: moran update rate for wild-type cells
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes.  
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should 
    match `length(μ)`.
- `ploidy::Int64 = 2`: cell ploidy 
"""
struct MoranInput <: SinglelevelInput
    N::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function MoranInput(;
    N = 1000,
    tmax = 15.0,
    moranrate = 1.0,
    moranincludeself = true,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)   
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return MoranInput(
        N, tmax, moranrate, moranincludeself, clonalmutations, μ, mutationdist, ploidy
    )
end

"""
    BranchingMoranInput <: SinglelevelInput

Input type for a single level simulation that grows by a branching process to `Nmax` cells and
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
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes. 
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should 
    match `length(μ)`.
- `ploidy::Int64 = 2`: cell ploidy
"""
struct BranchingMoranInput <: SinglelevelInput
    Nmax::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    birthrate::Float64
    deathrate::Float64
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function BranchingMoranInput(;
    Nmax = 1000,
    tmax = 15.0,
    moranrate = 1.0,
    moranincludeself = true,
    birthrate = moranrate,
    deathrate = 0.0,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)   
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return BranchingMoranInput(
        Nmax, tmax, moranrate, moranincludeself, birthrate, deathrate, clonalmutations, μ, 
            mutationdist, ploidy
    )
end
#endregion

#region Multi-level simulation inputs
abstract type MultilevelInput <: SimulationInput end

"""
    MultilevelBranchingInput <: MultilevelInput

Input type for a multilevel branching simulation that starts with a single cell in a single module. 
    
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
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes. 
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should 
    match `length(μ)`.
- `ploidy::Int64 = 2`
"""
struct MultilevelBranchingInput <: MultilevelInput
    modulesize::Int64
    maxmodules::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    asymmetricrate::Float64
    birthrate::Float64
    deathrate::Float64
    branchrate::Float64
    branchinitsize::Int64
    modulebranching::Symbol
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function MultilevelBranchingInput(;
    modulesize = 20,
    maxmodules = 1000,
    tmax = Inf,
    moranrate = 1.0,
    moranincludeself = true,
    asymmetricrate = 0.0,
    birthrate = maximum((moranrate, asymmetricrate)),
    deathrate = 0.0,
    branchrate  = 5.0,
    branchinitsize = 1,
    modulebranching = :split,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)   
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return MultilevelBranchingInput(
        modulesize, maxmodules, tmax, moranrate, moranincludeself, asymmetricrate, 
            birthrate, deathrate, branchrate, branchinitsize, modulebranching,
            clonalmutations, μ, mutationdist, ploidy
    )
end

"""
    MultilevelMoranInput <: MultilevelInput

Input type for a multilevel branching simulation that starts with `maxmodules` modules, each with
a single cell. 
    
Within module dynamics follows a branching process until `modulesize` is reached and then
switches to a Moran process. Module level dynamics follows a Moran process at rate 
`branchrate`.

# Keyword arguments:
- `modulesize::Int64 = 200`: maximum number of cells per module
- `maxmodules::Int64 = 10000`: maximum number of modules in population 
- `tmax::Float64 = Inf`: maximum time to run simulation
- `ploidy::Int64 = 2`: cell ploidy (per cell mutation rate is `ploidy * μ`)
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes. 
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells in branching phase
- `deathrate::Float64 = 0.0`: death rate for wild-type cells in branching phase
- `moranrate::Float64 = 1.0`: rate of Moran updating for wild-type cells in homeostasis
- `asymmetricrate::Float64 = 0.0`: rate of asymmetric updating for wild-type cells in homeostasis
- `branchrate::Float64 = 5.0`: rate at which homeostatic modules split to form new modules
- `branchinitsize::Int64 = 1`: number of cells sampled to form a 
    new module
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes. 
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should 
    match `length(μ)`.
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


struct MultilevelMoranInput <: MultilevelInput
    modulesize::Int64
    maxmodules::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    asymmetricrate::Float64
    birthrate::Float64
    deathrate::Float64
    branchrate::Float64
    branchinitsize::Int64
    modulebranching::Symbol
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function MultilevelMoranInput(;
    modulesize = 20,
    maxmodules = 1000,
    tmax = 15.0,
    moranrate = 1.0,
    moranincludeself = true,
    asymmetricrate = 0.0,
    birthrate = maximum((moranrate, asymmetricrate)),
    deathrate = 0.0,
    branchrate  = 5.0,
    branchinitsize = 1,
    modulebranching = :split,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)   
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return MultilevelMoranInput(
        modulesize, maxmodules, tmax, moranrate, moranincludeself, asymmetricrate, 
            birthrate, deathrate, branchrate, branchinitsize, modulebranching,
            clonalmutations, μ, mutationdist, ploidy
    )
end

"""
    MultilevelBranchingMoranInput <: MultilevelInput

Input type for a multilevel simulation of a homeostatic population that starts with a single cell 
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
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a 
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes. 
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new 
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should 
    match `length(μ)`.
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

struct MultilevelBranchingMoranInput <: MultilevelInput
    modulesize::Int64
    maxmodules::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    asymmetricrate::Float64
    birthrate::Float64
    deathrate::Float64
    branchrate::Float64
    branchinitsize::Int64
    modulebranching::Symbol
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function MultilevelBranchingMoranInput(;
    modulesize = 20,
    maxmodules = 1000,
    tmax = 15.0,
    moranrate = 1.0,
    moranincludeself = true,
    asymmetricrate = 0.0,
    birthrate = maximum((moranrate, asymmetricrate)),
    deathrate = 0.0,
    branchrate  = 5.0,
    branchinitsize = 1,
    modulebranching = :split,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)   
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return MultilevelBranchingMoranInput(
        modulesize, maxmodules, tmax, moranrate, moranincludeself, asymmetricrate, 
            birthrate, deathrate, branchrate, branchinitsize, modulebranching,
            clonalmutations, μ, mutationdist, ploidy
    )
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


function Base.show(io::IO, input::BranchingInput)
    @printf(io, "Single level branching process:\n")
    @printf(io, "    Maximum cells = %d\n", input.Nmax)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Birth rate = %.3f, death rate = %.3f\n", input.birthrate, input.deathrate)
    for i in 1:length(input.μ)
        @printf(io, "    %s: μ = %.3f\n", mutationdist_string(input.mutationdist[i]), input.μ[i])
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)
end 
  
function Base.show(io::IO, input::MoranInput)
    @printf(io, "Single level Moran process:\n")
    @printf(io, "    Maximum cells = %d\n", input.N)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Moran rate = %.3f \n", input.moranrate)
    for i in 1:length(input.μ)
        @printf(io, "    %s: μ = %.3f\n", mutationdist_string(input.mutationdist[i]), input.μ[i])
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)

end

function Base.show(io::IO, input::BranchingMoranInput)
    @printf(io, "Single level Branching -> Moran process:\n")
    @printf(io, "    Maximum cells = %d\n", input.Nmax)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Birth rate = %.3f, death rate = %.3f\n", input.birthrate, input.deathrate)
    @printf(io, "    Moran rate = %.3f\n", input.moranrate)
    for i in 1:length(input.μ)
        @printf(io, "    %s: μ = %.3f\n", mutationdist_string(input.mutationdist[i]), input.μ[i])
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)

end

function Base.show(io::IO, input::MultilevelInput)
    @printf(io, "Multilevel branching process:\n")
    @printf(io, "    Maximum modules = %d\n", input.maxmodules)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Module formation rate = %.2f\n", input.branchrate)
    @printf(io, "    Module formation mechanism = %s\n", input.modulebranching)
    @printf(io, "    Module size = %d\n", getmaxmodulesize(input))
    @printf(io, "    Number module founder cells = %d\n", input.branchinitsize)
    @printf(io, "    Birth rate = %.3f, death rate = %.3f\n", input.birthrate, input.deathrate)
    includestring = input.moranincludeself ? "include" : "exclude"
    @printf(io, "    Moran rate = %.3f (%s self)\n", input.moranrate, includestring)
    @printf(io, "    Asymmetric rate = %.3f\n", input.asymmetricrate)
    for i in 1:length(input.μ)
        @printf(io, "    %s: μ = %.3f\n", mutationdist_string(input.mutationdist[i]), input.μ[i])
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)
end

function mutationdist_string(mutationdist)
    if mutationdist == :fixed
        return "Fixed mutations"
    elseif mutationdist == :fixedtimedep
        return "Time-dependent fixed mutations"
    elseif mutationdist == :poisson
        return "Poisson distributed mutations"
    elseif mutationdist == :poissontimedep
        return "Time-dependent Poisson distibuted mutations"
    elseif mutationdist == :geometric
        return "Geometric distributed mutations"
    end
end