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

function sample_sims(multsim::MultiSimulation, n, rng::AbstractRNG = Random.GLOBAL_RNG)
    randindex = sample(rng, 1:length(multsim.output), n, replace = false)
    return [SomaticEvolution.get_simulation(multsim, i) for i in randindex]
end

function plot_example_VAF(multsim::MultiSimulation, n, rng::AbstractRNG=Random.GLOBAL_RNG; 
                            title = "", savedata = nothing)
    vafplots = []
    sampledsims = sample_sims(multsim, n, rng)
    for (i, simdata) in enumerate(sampledsims)
        VAFresult = getVAFresult(simdata, rng)
        df, fitcoef, rsq = fitinverse(VAFresult.sampledVAF, 0.12, 0.24)
        if savedata !== nothing
            CSV.write(savedata*"_$i.csv", df)
        end
        p1 = plotvaf(VAFresult, sampled=true, cumulative = false)
        p2 = plotinversevaf(VAFresult, sampled = true, fitcoef = fitcoef, rsq = rsq, 
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

function plot_example_VAF(VAFresult::VAFResult; title = "", savedata = nothing)

    df, fitcoef, rsq = fitinverse(VAFresult.sampledVAF, 0.12, 0.24)
    if savedata !== nothing
        CSV.write(savedata*".csv", df)
    end
    p1 = plotvaf(VAFresult, sampled=true, cumulative = false, savedata=savedata)
    p2 = plotinversevaf(VAFresult, sampled = true, fitcoef = fitcoef, rsq = rsq, 
        fmin=0.12, fmax=0.24, ylabel = "Cumulative number\nof mutations")
    label1 = latexstring("\\mu_\\textrm{pred} = $(round(Int,fitcoef))")
    label2 = latexstring("R^2 = $(round(rsq,digits=4))")
    annotate!([((0.1,0.96), text(label2,10,:left,:top))])
    annotate!([((0.1,0.85), text(label1,10,:left,:top))])
    p = plot(p1,p2,layout=(1,2),size = (600,300), plot_title=title, plot_fontsize=14,
            bottom_margin=10px,right_margin=10px, left_margin=10px,top_margin=15px)
    return p
      
end

function plot_singleplots(input1, input2, parentdir="", rng::AbstractRNG=Random.GLOBAL_RNG)

    dir = parentdir*"singleplots/"
    mkpath(dir)

    simdata1 = run1simulation(input1, rng)
    simdata2 = run1simulation(input2, rng)
    VAFresult1 = getVAFresult(simdata1, rng)
    VAFresult2 = getVAFresult(simdata2, rng)


    p1 = plot_example_VAF(VAFresult1, title = "No subclones, neutral tumour growth",
                            savedata=dir*"neutral")
    p2 = plot_example_VAF(VAFresult2, title = "One subclone with s = 1, t1 = 6",
                            savedata=dir*"oneclone")

    saveinput(dir*"neutral_input.txt", simdata1)
    saveinput(dir*"oneclone_input.txt", simdata2)

    savefig(p1, dir*"neutral_expVAF.pdf")
    savefig(p2, dir*"oneclone_expVAF.pdf")
    return (results=(simdata1, simdata2), plots=(p1,p2))
end

function plot_multiplots(input1, input2, parentdir="", rng::AbstractRNG=Random.GLOBAL_RNG)

    dir = parentdir*"multiplots/"
    mkpath(dir)

    results = multiplesimulations(100, input1, input2, rng = rng)

    p1, sims1 = plot_example_VAF(results[1], 3, rng, title = "No subclones, neutral tumour growth",
                        savedata=dir*"neutral")
    p2, sims2 = plot_example_VAF(results[2], 3, rng, title = "One subclone with s = 1, t1 = 6",
                        savedata=dir*"oneclone")
    df = getstats(results,:selection,:numclones)
    p = plot(plot_μ_error_vs_n(df),plot_r2_vs_n(df))
    savefig(p1, dir*"neutral_expVAF.pdf")
    savefig(p2, dir*"oneclone_expVAF.pdf")
    savefig(p, dir*"exp_stats.pdf")
    saveinput(dir*"neutral_input.txt", results[1], 1)
    for (i,sim) in enumerate(sims2)
        saveinput(dir*"oneclone_input_$i.txt", sim)
    end
    return (results=results, plots=(p1,p2,p))
end

rng = MersenneTwister(13)
μ = 100
b = log(2)

input1 = BranchingInput(Nmax=10000,numclones=0,μ=μ,clonalmutations=2*μ)
input2 = BranchingInput(Nmax=10000,numclones=1,μ=μ,clonalmutations=2*μ,
                selection = [1.0], tevent = [6.0])

sim = run1simulation(input2, rng)
results1, plots1 = plot_singleplots(input1, input2, "data/")
results2, plots2 = plot_multiplots(input1, input2, "data/");