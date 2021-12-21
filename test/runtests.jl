push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using SomaticEvolution
using Test

tests = ["init","multilevel"]

@testset "SomaticEvolution.jl" begin
    for test in tests
        @testset "$test" begin
            include(test*".jl")
        end
    end
end

