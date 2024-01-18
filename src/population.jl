abstract type AbstractPopulation end

struct Population{T<:AbstractModule} <: AbstractPopulation
    homeostatic_modules::Vector{T}
    growing_modules::Vector{T}
    subclones::Vector{Subclone}
end

struct SinglelevelPopulation{T<:AbstractModule} <: AbstractPopulation
    singlemodule::T
    subclones::Vector{Subclone}
end

Base.length(population::SinglelevelPopulation) = length(population.singlemodule)
Base.length(population::Population) = length(population.homeostatic_modules) + length(population.growing_modules)
Base.iterate(population::Population) = iterate(Base.Iterators.flatten((population.homeostatic_modules, population.growing_modules)))
Base.iterate(population::Population, state) = iterate(Base.Iterators.flatten((population.homeostatic_modules, population.growing_modules)), state)

function Base.getindex(population::Population, i)
    Nhom = length(population.homeostatic_modules)
    if i <= Nhom 
        return Base.getindex(population.homeostatic_modules, i)
    else
        return Base.getindex(population.growing_modules, i - Nhom)
    end
end

Base.getindex(population::Population, a::Vector{Int64}) = map(i -> getindex(population, i), a)

function Base.setindex(population::Population, v, i)
    Nhom = length(population.homeostatic_modules)
    if i <= Nhom 
        return Base.setindex(population.homeostatic_modules, v, i)
    else
        return Base.setindex(population.growing_modules, v, i - Nhom)
    end
end

function Base.firstindex(population::Population)
    if length(population.homeostatic_modules) != 0 
        return firstindex(population.homeostatic_modules)
    else
        return firstindex(population.growing_modules)
    end
end

function Base.lastindex(population::Population)
    if length(population.growing_modules) != 0 
        return lastindex(population.growing_modules)
    else 
        return lastindex(population.homeostatic_modules)
    end
end

function Population(homeostatic_modules::Vector{T}, growing_modules::Vector{T}, birthrate, deathrate, moranrate, asymmetricrate) where T
    return Population{T}(
        homeostatic_modules, 
        growing_modules,
    Subclone[Subclone(
        1,
        0,
        0.0,
        sum(length.(homeostatic_modules)) + sum(length.(growing_modules)),
        birthrate,
        deathrate,
        moranrate,
        asymmetricrate
    )]
    )
end

function SinglelevelPopulation(singlemodule::T, birthrate, deathrate, moranrate, asymmetricrate) where T
    return SinglelevelPopulation{T}(
        singlemodule,
    Subclone[Subclone(
        1,
        0,
        0.0,
        length(singlemodule),
        birthrate,
        deathrate,
        moranrate,
        asymmetricrate
    )]
    )
end

allmodules(population) = vcat(population.homeostatic_modules, population.growing_modules)

function number_cells_by_subclone(modules, nsubclones)
    ncells = zeros(Int64, nsubclones)
    for mod in modules
        for cell in mod.cells
            ncells[getclonetype(cell)] += 1
        end
    end
    return ncells
end

function number_cells_in_subclone(modules, subcloneid)
    return sum(getclonetype(cell) == subcloneid for mod in modules for cell in mod.cells)
end

function Base.show(io::IO, population::Population{T}) where T
    Base.show(io, Population{T})
    @printf(io, ": \n    %d growing modules", length(population.growing_modules))
    @printf(io, "\n    %d homeostatic modules", length(population.homeostatic_modules))
    @printf(io, "\n    %d subclones", length(filter(x -> x.size > 0, population.subclones)))
end

function Base.show(io::IO, population::SinglelevelPopulation{T}) where T
    Base.show(io, SinglelevelPopulation{T})
    mod = population.singlemodule
    @printf(io, ": \n    %d cells ", length(mod))
    printmodule(io, mod)
    @printf(io, " (t = %.2f)", age(mod))
end

function Base.show(io::IO, mod::T) where T<:AbstractModule
    Base.show(io, T)
    @printf(io, ": \n    %d cells ", length(mod))
    printmodule(io, mod)
    @printf(io, " (t = %.2f)", age(mod))


end

function printmodule(io::IO, mod::AbstractModule; maxsubclone=nothing)
    if length(mod) == 0
        @printf(io, "[]")
    else
        subclonelist = sort!([getclonetype(cell) for cell in mod])
        unique_subclones = unique(subclonelist)
        @printf(io, "[")
        if length(unique_subclones) == 1 
            @printf(io, "%d (%d)", length(subclonelist), unique_subclones[1])
        else
            for subcloneid in unique_subclones[1:end-1]
                @printf(io, "%d (%d), ", count(x->x==subcloneid, subclonelist), subcloneid)
            end
            subcloneid = unique_subclones[end]
            @printf(io, "%d (%d)", count(x->x==subcloneid, subclonelist), subcloneid)
        end
        @printf(io, "]")
    end
end

# function printmodule(io::IO, mod::AbstractModule; maxsubclone=nothing)
#     if length(mod) == 0
#         @printf(io, "[]")
#     else
#         subclonelist = [getclonetype(cell) for cell in mod]
#         maxsubclone = isnothing(maxsubclone) ? maximum(subclonelist) : maxsubclone
#         subclone_counts = counts(subclonelist, 1:maxsubclone)
#         @printf(io, "[%d", subclone_counts[1])
#         if length(subclone_counts) > 1
#             for c in subclone_counts[2:end]
#                 @printf(io, ", %d", c)
#             end
#         end
#         @printf(io, "]")
#     end
# end
