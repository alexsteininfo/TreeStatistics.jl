using SomaticEvolution
using Test
using Random
using StatsBase

tests = ["initialisation","multilevel","simulations"]

@testset "SomaticEvolution.jl" begin
    for test in tests
        @testset "$test" begin
            include(test*".jl")
        end
    end
end

