using Revise
using SomaticEvolution
using Test
using Random
using StatsBase
using AbstractTrees
using Distributions

tests = [
    "multilevel_tree",
    "initialisation",
    "multilevel",
    "simulations",
    "process_mutations",
    "treesimulations",
    "multilevel_selection"
]

@testset "SomaticEvolution.jl" begin
    for test in tests
        @testset "$test" begin
            include(test*".jl")
        end
    end
end