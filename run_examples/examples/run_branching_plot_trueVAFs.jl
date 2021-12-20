push!(LOAD_PATH, "/Users/renton02/reps/somatic-evolution/SomaticEvolution")

using Revise
using Random
using StatsBase
using SomaticEvolution
using Plots
using Plots.PlotMeasures
using StatsPlots
using LaTeXStrings
using CSV

# Plots.scalefontsizes(1/1.2)

function plotsim(simdata, dir)
    mkpath(dir)
    _, fitcoef1, rsq1 = fitinverse(simdata.sampled.VAF, 0.12, 0.24)
    _, fitcoef2, rsq2 = fitinverse(simdata.output.trueVAF, 0.01, 0.12)
    p1a = plotvaf(simdata, sampled=true, cumulative = false)
    p1b = plotinversevaf(simdata, sampled = true, fitcoef = fitcoef1, rsq = rsq1, 
            fmin=0.12, fmax=0.24, ylabel = "Cumulative number\nof mutations", 
            leftmargin=20px, rightmargin=12px, lw=2)
    label1 = latexstring("\\mu_\\textrm{pred} = $(round(Int,fitcoef1[1]))")
    label2 = latexstring("R^2 = $(round(rsq1,digits=4))")
    annotate!([((0.1,0.96), text(label2,16,:left,:top))])
    annotate!([((0.1,0.84), text(label1,16,:left,:top))])
    # p1 = plot(p1a,p1b,layout=(1,2), size = (600,320))

    p2a = plotvaf(simdata, sampled=false, cumulative = false)
    p2b = plotinversevaf(simdata, sampled = false, fitcoef = fitcoef2, rsq = rsq2, 
            fmin=0.01, fmax=0.12, ylabel = "Cumulative number\nof mutations", 
            leftmargin=20px, rightmargin=12px, lw=2)
    label1 = latexstring("\\mu_\\textrm{pred} = $(round(Int,fitcoef2[1]))")
    label2 = latexstring("R^2 = $(round(rsq2,digits=4))")
    annotate!([((0.1,0.96), text(label2,16,:left,:top))])
    annotate!([((0.1,0.84), text(label1,16,:left,:top))])
    # p2 = plot(p2a,p2b,layout=(1,2), size=(600, 320))


    savefig(p1a, dir*"branching_sampled_VAF.pdf")
    savefig(p1b, dir*"branching_sampled_invVAF.pdf")
    savefig(p2a, dir*"branching_true_VAF.pdf")
    savefig(p2b, dir*"branching_true_invVAF.pdf")
end

rng = MersenneTwister(87291)
input1 = InputParameters{BranchingInput}(Nmax=10000,numclones=0,μ=μ,clonalmutations=2*μ)
simdata1 = run1simulation(input1, rng)

rng = MersenneTwister(81)
input2 = InputParameters{BranchingInput}(Nmax=10000,numclones=1,μ=μ,clonalmutations=2*μ,
            selection=[1], tevent=[6])
simdata2 = run1simulation(input2, rng)

dir = "presplots/"
plotsim(simdata1, dir*"neutral/")
plotsim(simdata2, dir*"oneclone/")