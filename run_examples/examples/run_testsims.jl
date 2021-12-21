push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots

rng = MersenneTwister(14)

IP1 = InputParameters{BranchingMoranInput}(Nmax=100, numclones=0, μ=100, fixedmu=false, 
    bdrate=log(2), clonalmutations=200, tmax=20)

IP2 = InputParameters{BranchingInput}(Nmax=100, numclones=0, μ=100, fixedmu=false, b=log(2),
    clonalmutations=200)

IP3 = InputParameters{MoranInput}(N=100, numclones=0, μ=100, fixedmu=false, bdrate=log(2),
    clonalmutations=200, tmax=20)

IP4 = InputParameters{BranchingMoranInput}(Nmax=100, tmax=20, numclones=2, μ=100, fixedmu=false, 
    bdrate=log(2), clonalmutations=200, selection=[1,2], tevent=[5,6])

IP5 = InputParameters{BranchingInput}(Nmax=100,numclones=2,μ=100,clonalmutations=200,
    selection=[1,2], tevent=[4,5])

IP6 = InputParameters{MoranInput}(N=100,tmax=20,numclones=2,μ=100,clonalmutations=200,
    selection=[1,2], tevent=[5,6])

plots = []
for input in (IP1, IP2, IP3, IP4, IP5, IP6)
    sresult = run1simulation(input, rng)
    println(sresult)
    push!(plots, plotpopulation(sresult))
end

display(plot(plots..., layout=(3,2), link=:all, legend=false))
