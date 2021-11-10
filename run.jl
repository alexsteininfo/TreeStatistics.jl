push!(LOAD_PATH, "/Users/jessie/git_reps/somatic-evolution/SomaticEvolution")

using Revise
using Random
using SomaticEvolution
using Plots

rng = MersenneTwister(13)

input = InputParameters{BranchingInput}(Nmax=10000,numclones=0,μ=200,clonalmutations=400)
simdata = run1simulation(input,rng)

@userplot struct PlotVAF{T<:Tuple{Simulation}}
    args::T
end

@recipe function f(pv::PlotVAF; sampled = true, cumulative = false)
    df = sampled ? pv.args[1].sampled.df : gethist(pv.args[1].output.trueVAF)
    VAF = df[!,:VAF]
    freq = cumulative ? df[!,:cumfreq] : df[!,:freq]

    @series begin
        seriestype --> :bar
        linecolor --> :darkslategrey
        fillcolor --> :darkslategrey
        markerstrokecolor --> :white
        VAF, freq
    end

    if length(pv.args[1].output.clonefreq) > 0
        xint = pv.args[1].output.clonefreq./2 * pv.args[1].input.cellularity

        @series begin
            seriestype --> :vline
            legend --> false
            fillcolor --> :darkred
            linewidth --> 3
            xint
        end
    end
    yaxis --> cumulative ? "Cumulative number of mutations" : "Number of mutations"
    xaxis --> "VAF"
    legend --> false
    grid --> false
    ()
end

@userplot struct PlotInverseVAF{T<:Tuple{Simulation}}
    args::T
end


@recipe function f(piv::PlotInverseVAF; sampled = true, cumulative = true)
    df = sampled ? piv.args[1].sampled.df : gethist(piv.args[1].output.trueVAF)
    VAF = df[!,:VAF]
    freq = cumulative ? df[!,:cumfreq] : df[!,:freq]

    @series begin
        seriestype --> :scatter
        markercolor --> :darkslategrey
        markerstrokecolor --> :white
        1 ./ VAF, freq
    end
    μ = piv.args[1].input.siminput.μ 
    β = 1 - piv.args[1].input.siminput.d/piv.args[1].input.siminput.b
    ploidy = piv.args[1].input.ploidy
    @series begin
        seriestype --> :line
        linecolor --> :red
        1 ./ VAF, μ/β .* (1 ./ VAF .- ploidy)
    end

    if length(piv.args[1].output.clonefreq) > 0
        xint = piv.args[1].output.clonefreq./2 * piv.args[1].input.cellularity


        @series begin
            seriestype --> :vline
            fillcolor --> :darkred
            linewidth --> 3
            xint
        end
    end
    yaxis --> cumulative ? "Cumulative number of mutations" : "Number of mutations"
    xaxis --> "Inverse VAF"
    xmin,xmax = 1/0.24,1/0.12
    xlims --> (xmin,xmax)
    ylims --> (0,2000)
    xticks --> ([1/0.24, 1/0.16, 1/0.12],["1/0.24", "1/0.16", "1/0.12"])
    legend --> false
    grid --> false
    ()
end

p1 = plotvaf(simdata,legend=false, grid=false, sampled=true, cumulative = true)
p2 = plotinversevaf(simdata, sampled = true)
plot(p1,p2,layout=(2,1))