push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using SomaticEvolution
using Test

tests = []

@testset "SomaticEvolution.jl" begin
    for test in tests
        include(test*".jl")
    end
end

