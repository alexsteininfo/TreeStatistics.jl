abstract type AbstractPopulation end

struct Population{T<:AbstractModule} <: AbstractPopulation
    homeostatic_modules = Vector{T}
    growing_modules = Vector{T}
    subclones = Vector{Subclone}
end

Base.length(population::Population) = length(population.homeostatic_modules) + length(population.growing_modules)
Base.iterate(population::Population) = iterate(vcat(population.homeostatic_modules, population.growing_modules))
Base.iterate(population::Population, state) = iterate(vcat(population.homeostatic_modules, population.growing_modules), state)

function Population(growing_modules::Vector{T}, birthrate, deathrate, moranrate, asymmetricrate) where T
    return Population{T}(
        T[], 
        growing_modules,
    Subclone[Subclone(
        1,
        0.0,
        sum(length.(growing_modules)),
        birthrate,
        deathrate,
        moranrate,
        asymmetricrate
    )]
    )
end

function Population(homeostatic_modules::Vector{T}, growing_modules::Vector{T}, birthrate, deathrate, moranrate, asymmetricrate) where T
    return Population{T}(
        homeostatic_modules, 
        growing_modules,
    Subclone[Subclone(
        1,
        0.0,
        sum(length.(growing_modules)) + sum(length.(growing_modules)),
        birthrate,
        deathrate,
        moranrate,
        asymmetricrate
    )]
    )
end

allmodules(population) = vcat(population.homeostatic_modules, population.growing_modules)

function getwildtyperates(population::Population)
    return (
        birthrate = population.subclones[1].birthrate,
        deathrate = population.subclones[1].deathrate,
        moranrate = population.subclones[1].moranrate,
        asymmetricrate = population.subclones[1].asymmetricrate
    )
end

