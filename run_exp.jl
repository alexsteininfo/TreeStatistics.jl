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

rng = MersenneTwister()
μ = 100
b = (log(2))

input1 = InputParameters{BranchingInput}(Nmax=10000,numclones=0,μ=μ,clonalmutations=2*μ)
input2 = InputParameters{BranchingInput}(Nmax=10000,numclones=1,μ=μ,clonalmutations=2*μ,
                    selection = [1.0], tevent = [6.0])


results = multiplesimulations(100, input1, input2, rng = rng)
# p1, sims1 = plot_example_VAF(results[1], 3, rng, title = "No subclones, neutral tumour growth")
# p2, sims2 = plot_example_VAF(results[2], 3, rng, title = "One subclone with s = 1, t1 = 6")

# df = getstats(results,:selection,:numclones)
# p = plot(plot_μ_error_vs_n(df),plot_r2_vs_n(df))

# savefig(p1, "neutral_expVAF.pdf")
# savefig(p2, "oneclone_expVAF.pdf")
# savefig(p, "exp_stats.pdf")
