push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots

rng = MersenneTwister(12)

IP = InputParameters{MultilevelInput}(Nmax=100, numclones=0, μ=100, fixedmu=false, bdrate=log(2),
    clonalmutations=0, pop_age=10, branchrate=1/5, branchinitsize=10)

pop = multilevel_simulation(IP, rng=rng)

