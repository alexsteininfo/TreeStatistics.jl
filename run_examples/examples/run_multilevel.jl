push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots; gr()


# parameters: 
# bdrate = 1, in homeostatic population assume birth rate of 1 per cell per day
# b = 10, d = 0, in growing population assume birth rate is 10x faster and no deaths
# clonalmutations = 0, initial cell has no mutations
# pop_age = 17*365, total time (in days) for population growth
# branchrate = 3/365 - 8/365, 3-8 per year
# modulesize = 20-200 cells, homeostatic size for each module
# branchfraction = 0.05-0.5, proportion of cells to sample from module to form new module
# μ = (1e-8 - 1e-9) x 260e6 x 2, mutation rate per base pair x genome size x ploidy

function main(modulenumber, age, seed, showprogress=false)
    rng = MersenneTwister(seed)
    input = MultilevelInput(modulesize=20, numclones=0, fixedmu=false, b=1, d=0,
        bdrate=1, clonalmutations=0, pop_age=age*365, branchrate=3/365, branchfraction=0.2, 
        μ=1e-8*260e6*2)
    pop = multilevel_simulation(input, rng=rng, maxmodules=modulenumber, showprogress=showprogress)
    return pop
end

function main_fast(modulenumber, age, seed, showprogress=false)
    rng = MersenneTwister(seed)
    input = MultilevelInput(modulesize=20, numclones=0, fixedmu=false, b=1, d=0,
        bdrate=1, clonalmutations=0, pop_age=age*365, branchrate=3/365, branchfraction=0.2, 
        μ=1e-8*260e6*2)
    pop = SomaticEvolution.multilevel_simulation_fast(input, rng=rng, maxmodules=modulenumber, showprogress=showprogress)
    return pop
end

function plotpairwise(pop, idx=nothing, n=16; rng=Random.GLOBAL_RNG)
    if isnothing(idx)
        idx = sample(rng, 1:length(pop), n, replace=false)
    end
    pfd = pairwise_fixed_differences(pop, idx)
    return pairwiseplot(pfd), idx
end

function plotVAFs(pop, idx=nothing, n=16; rng=Random.GLOBAL_RNG)
    if isnothing(idx)
        idx = sample(rng, 1:length(pop), n, replace=false)
    end
    pvec = [plotvaf(get_simulation(pop,i), xlabel="", ylabel="") for i in idx]
    p2 = plot(pvec..., layout=(4,4), legend=false)
end