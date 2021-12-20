push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots

rng = MersenneTwister(14)

input1 = InputParameters{BranchingMoranInput}(Nmax=100, numclones=0, μ=100, fixedmu=false, 
    bdrate=log(2), clonalmutations=200, tmax=10)

input2 = InputParameters{BranchingInput}(Nmax=100, numclones=0, μ=100, fixedmu=false, b=log(2),
    clonalmutations=200)

input3 = InputParameters{MoranInput}(N=100, numclones=0, μ=100, fixedmu=false, bdrate=log(2),
    clonalmutations=200, tmax=10)

# p1 = plotvaf(simdata, sampled=true, cumulative=false)
# df, fitcoef, rsq = fitinverse(simdata.output.trueVAF, 0.0, 0.05, 0.005, false)
# p2 = plotinversevaf(simdata, sampled=false, cumulative=false, fitcoef=fitcoef,
#                     fmin=0, fmax=0.12, fstep=0.005, dataseries=:scatter)

plots = []
for input in (input1, input2, input3)
    sresult = run1simulation(input, rng)
    push!(plots, plotpopulation(sresult))
end

display(plot(plots..., layout=(3,1), link=:all))
