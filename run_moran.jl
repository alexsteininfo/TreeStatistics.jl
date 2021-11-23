push!(LOAD_PATH, "/Users/jessie/git_reps/somatic-evolution/SomaticEvolution")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots
using StatsPlots
using LaTeXStrings

function plot_μ_error_vs_n(df)
    # p = @df df violin(:numclones, :μ_error, linewidth = 0, barwidth = 0.5)
    p = @df df boxplot(:numclones, :μ_error, bar_width = 0.3, legend = false)
    xticks!([0,1,2])
    xlabel!("number of subclones")
    ylabel!("% error on μ estimate")
    return p
end

function plot_r2_vs_n(df)
    # p = @df df violin(:numclones, :r2, linewidth = 0, barwidth = 0.5)
    p = @df df boxplot(:numclones, :r2, bar_width = 0.3, legend = false)
    xticks!([0,1,2])
    xlabel!("number of subclones")
    ylabel!("R^2")
    return p
end

function sample_sims(multsim::MultiSimulation, n, rng::AbstractRNG = MersenneTwister())
    simdatalist = Simulation[]
    randindex = sample(rng, 1:length(multsim.output), n, replace = false)
    for i in randindex
        simdata = Simulation(multsim.input, multsim.output[i], multsim.sampled[i])
        push!(simdatalist, simdata)
    end
    return simdatalist
end

function plot_example_VAF(multsim, n, rng::AbstractRNG = MersenneTwister(); title = "")
    vafplots = []
    sampledsims = sample_sims(multsim, n, rng)
    for (i, simdata) in enumerate(sampledsims)
        df, fitcoef, rsq = fitinverse(simdata.sampled.VAF, 0.12, 0.24)
        p1 = plotvaf(simdata, sampled=true, cumulative = false)
        p2 = plotinversevaf(simdata, sampled = true, fitcoef = fitcoef, rsq = rsq, 
            fmin=0.12, fmax=0.24)
        label1 = latexstring("\\mu_\\textrm{pred} = $(round(Int,fitcoef))")
        label2 = latexstring("R^2 = $(round(rsq,digits=4))")
        annotate!([((0.1,0.96), text(label2,10,:left,:top))])
        annotate!([((0.1,0.7), text(label1,10,:left,:top))])
        p = plot(p1,p2,layout=(1,2))
        if i != (n+(n%2))/2 ylabel!(" ") end
        push!(vafplots,p)
    end
    p = plot(vafplots..., layout = (n,1), link=:x, plot_title = title, 
        plot_titlefontsize = 14, size = (550,500))
    return p, sampledsims
      
end

rng = MersenneTwister(13)
μ = 100
b = (log(2))

input1 = InputParameters{MoranInput}(N=100,numclones=0,μ=μ,clonalmutations=100, tmax=500)
# input2 = InputParameters{MoranInput}(N=100,numclones=1,μ=μ,clonalmutations=2*μ, tmax=5,
#                     selection = [10.0], tevent = [4.0])
simdata1 = run1simulation(input1, rng)
plotvaf(simdata1, sampled=false, cumulative=true)
# simdata2 = run1simulation(input2, rng)

# simtracker = SomaticEvolution.moranprocess(10, log(2), 5, 1, rng, clonalmutations=0,
#     numclones=1, selection=[10.0], tevent=[4.0])

# simtracker, _ , simresults = SomaticEvolution.processresults!(simtracker,100,1,μ,false,2*μ,2,0,1,rng)