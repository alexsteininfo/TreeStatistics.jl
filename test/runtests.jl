push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using SomaticEvolution
using Test

moran_neutralIP = InputParameters{MoranInput}(N=100, numclones=0, μ=100, fixedmu=false, clonalmutations=0, tmax=20/log(2))
inputs = [moran_neutralIP]

@testset "SomaticEvolution.jl" begin
    for IP in inputs
        simtracker = SomaticEvolution.initializesim(IP.siminput)
        @test length(simtracker.subclones) == IP.siminput.numclones
    end
end

