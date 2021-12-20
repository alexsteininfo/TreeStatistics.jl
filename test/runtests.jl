push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using SomaticEvolution
using Test

tests = ["moran_tests.jl"]

@testset "SomaticEvolution.jl" begin
    for test in tests
        include(test)
    end
end

