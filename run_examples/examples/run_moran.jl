push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/src")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots
using StatsPlots
using LaTeXStrings


rng = MersenneTwister(11)
μ = 10
b = (log(2))

# input1 = InputParameters{BranchingInput}(Nmax=10000,numclones=0,μ=μ,clonalmutations=0)
# simdata1 = run1simulation(input1, rng)

input2 = InputParameters{MoranInput}(N=10000, numclones=0, μ=100, fixedmu=true, clonalmutations=0, tmax=1000/log(2))
simdata2 = run1simulation(input2, rng)

p1 = plotvaf(simdata2, sampled=false, cumulative=false, fstep=0.01)

p2 = plotinversevaf(simdata2, sampled=false, cumulative=false, 
                    fmin=0.0001, fmax=0.12, fstep=0.0001, dataseries=:scatter)


for p in (p1,p2)
    display(p)
end

# simdata2 = run1simulation(input2, rng)

# simtracker = SomaticEvolution.moranprocess(10, log(2), 5, 1, rng, clonalmutations=0,
#     numclones=1, selection=[10.0], tevent=[4.0])

# simtracker, _ , simresults = SomaticEvolution.processresults!(simtracker,100,1,μ,false,2*μ,2,0,1,rng)