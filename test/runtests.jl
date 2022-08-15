using SomaticEvolution
using Test
using Random
using StatsBase
using AbstractTrees

tests = [
    "multilevel_tree",
    "initialisation",
    "multilevel",
    "simulations",
    "testio",
    "process_mutations",
    "analysis",
    "treesimulations"
]

@testset "SomaticEvolution.jl" begin
    for test in tests
        @testset "$test" begin
            include(test*".jl")
        end
    end
end

