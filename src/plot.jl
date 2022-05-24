""" plotvaf(VAF, [subclonefreqs=[]]; <keyword arguments>)
"""
# plotvaf

@userplot PlotVAF

@recipe function f(pv::PlotVAF; cumulative = false, fstep = 0.01, sampled=true)
    VAFresult = pv.args[1]
    VAF = sampled ? VAFresult.sampledVAF : VAFresult.trueVAF
    df = gethist(VAF, fstep = fstep)
    VAF = (df[!,:VAF] .*2 .- fstep) ./ 2 #set x values to middle of bins
    freq = cumulative ? df[:,:cumfreq] : df[:,:freq]

    if get(plotattributes, :yscale, nothing) == :log10
        VAF, freq = VAF[freq .> 0], freq[freq .>0]
    end
    yguide = cumulative ? "Cumulative number of mutations" : "Number of mutations"

    @series begin
        seriestype --> :bar
        linecolor --> :white
        linewidth --> 0.2
        fillcolor --> :darkslategrey
        markerstrokecolor --> :white
        VAF, freq
    end

    if length(VAFresult.subclonefreq) > 0
        subclonefreq = VAFresult.subclonefreq ./VAFresult.input.ploidy * VAFresult.cellularity
    
        @series begin
            seriestype --> :vline
            legend --> false
            fillcolor --> :darkred
            linewidth --> 2
            legend --> false
            grid --> false
            subclonefreq
        end
    end
    yguide --> yguide
    xguide --> "VAF"
    legend --> false
    grid --> false
    ()
end

""" plotinversevaf(VAF; <keyword arguments>)
"""
plotvaf

@userplot PlotInverseVAF

@recipe function f(pv::PlotInverseVAF; fmin = 0.12, fmax = 0.24, fit_fmax=fmax, sampled=true,
                    fstep = 0.001, fitcoef = nothing, cumulative = true, dataseries = :line)

    VAFresult = pv.args[1]
    VAF = sampled ? VAFresult.sampledVAF : VAFresult.trueVAF
    df = gethist(VAF, fmin = fmin, fmax = fmax, fstep = fstep) 
    VAF = df[!,:VAF]
    freq = cumulative ? df[!,:cumfreq] : df[!,:freq]
    yguide = cumulative ? "Cumulative numberof mutations" : "Number of mutations"

    @series begin
        yguide --> yguide
        xguide --> "Inverse VAF"
        seriestype --> dataseries
        1 ./ VAF, freq
    end

    if fitcoef !== nothing
        if length(fitcoef) == 1
            m = fitcoef[1]
            c = - m/fit_fmax
        else
            c, m = fitcoef
        end
        x, y = 1 ./ VAF, m ./ VAF .+ c
        x, y = x[y .>= 0], y[y .>= 0]
        @series begin
            yguide --> yguide
            xguide --> "Inverse VAF"
            seriestype --> :line
            linecolor --> :red
            linestyle --> :dash
            x, y
        end
    end

    yguide --> yguide
    xguide --> "Inverse VAF"
    fvals = [df[end, :VAF], 2/(1/df[1, :VAF]+ 1/df[end, :VAF]), df[1, :VAF]] 
    xticks --> (1 ./ fvals, map(x -> "1/" * x, string.(round.(fvals, digits = 2))))
    legend --> false
    grid --> false
    ()
end



@recipe function f(output::ModuleTracker)
    @series begin
        yguide --> "Population size"
        xguide --> "Time"
        seriestype --> :line
        output.tvec, output.Nvec
    end
end

@recipe function f(simulation::Simulation; fitexponential=false)
    yguide --> "Number of cells"
    xguide --> "Time"
    @series begin
        seriestype --> :line
        legend --> false
        label --> "simulation"
        grid --> false
        simulation.output.tvec, simulation.output.Nvec
    end
    if fitexponential
        df = DataFrame(
            :popsizelog => log.(simulation.output.Nvec), 
            :time => simulation.output.tvec
        )
        lmfit = fit(LinearModel, @formula(popsizelog ~ time + 0), df)
        m = coef(lmfit)[1]
        @series begin
            seriestype --> :line
            legend --> true
            legend_position --> :topleft
            label --> "b = $(round(m, digits=5))"
            simulation.output.tvec, exp.(m*simulation.output.tvec)
        end
    end
end

@recipe function f(multisim::MultiSimulation; plottype=:popsize, fitexponential=false, tstep=nothing)
    if plottype == :modulesize
        yguide --> "Module size"
        xguide --> "Time"
        for moduletracker in multisim
            @series begin
                seriestype --> :line
                legend --> false
                moduletracker.tvec, moduletracker.Nvec
            end
        end
    elseif plottype == :popsize
        time, popsize = newmoduletimes(multisim), 1:length(multisim)
        yguide --> "Number of modules"
        xguide --> "Time"
        @series begin
            seriestype --> :line
            legend --> false
            label --> "simulation"
            time, popsize
        end
        if fitexponential
            df = DataFrame(
                :popsizelog => log.(popsize), 
                :time => time
            )
            lmfit = fit(LinearModel, @formula(popsizelog ~ time + 0), df)
            m = coef(lmfit)[1]
            @series begin
                seriestype --> :line
                legend --> true
                legend_position --> :topleft
                label --> "r = $(round(m, digits=5))"
                time, exp.(m*time)
            end
        end

    elseif plottype == :cellpopsize
        tstep = isnothing(tstep) ? 1 / multisim.input.bdrate : tstep
        yguide --> "Number of cells"
        xguide --> "Time"
        @series begin
            seriestype --> :line
            legend --> false
            cellpopulationsize(multisim, tstep)
        end
    end
    
end

@userplot PairwisePlot

@recipe function f(pp::PairwisePlot; sampleids=nothing)
    if length(pp.args) == 1
        z = pp.args[1]
        n = isnothing(sampleids) ? size(z)[1] : length(sampleids)
        x, y = 1:n, 1:n
    elseif length(pp.args) == 3
        x, y, z = pp.args
    end
    if sampleids !== nothing
        z = z[sampleids, sampleids] 
    end
    yguide --> "Module"
    xguide --> "Module"
    title --> "Pairwise fixed differences"
    @series begin
        grid --> false
        seriestype := :heatmap
        seriescolor --> cgrad(:bone, rev=true)
        xticks --> x
        yticks --> y
        aspectratio --> 1
        size --> (500,510)
        x, y, z
    end
end

@userplot ModuleSharedPlot

@recipe function f(msp::ModuleSharedPlot)
    @series begin
        grid --> false
        seriestype := :bar
        xguide --> "Number of modules"
        yguide --> "Number fixed mutations"
        legend --> false
        msp.args[1]
    end
end

@userplot PairwiseDistributionPlot

@recipe function f(pdp::PairwiseDistributionPlot; sampleids=nothing, samplesize=nothing,
    showclonal=false, xticks=:auto)

    if typeof(pdp.args[1]) <: Dict
        pfddata = reduce(vcat, [fill(key, val) for (key,val) in pdp.args[1]])
    else
        pfdmatrix = pdp.args[1]
        n = size(pfdmatrix)[1]
        if sampleids === nothing && samplesize !== nothing
            sampleids = sample(1:n, samplesize, replace=false)
        end
        if sampleids !== nothing
            pfdmatrix = pfdmatrix[sampleids, sampleids]
            n = length(sampleids)
        end
        pfddata = convert_pfdmatrix_to_vector(pfdmatrix)
    end
    @series begin
        grid --> false
        seriestype := :stephist
        xguide --> "Number of pairwise fixed differences"
        yguide --> "Frequency"
        fillcolor --> :lightgrey
        fillrange --> 0
        linecolor --> :black
        legend --> false
        titlefontsize -->10
        titlelocation --> :left
        normalize --> :pdf
        pfddata
    end
    if showclonal && length(pdp.args) > 1
        @series begin
            legend --> false
            grid --> false
            xguide --> "Number of pairwise fixed differences"
            yguide --> "Frequency"
            titlefontsize -->10
            titlelocation --> :left
            seriestype --> :vline
            fillcolor --> :red
            linewidth --> 2
            [mean(pdp.args[2])]
        end
    end
end

function convert_pfdmatrix_to_vector(pfdmatrix)
    n = size(pfdmatrix)[1]
    vals = Int64[]
    for i in 1:n
        for j in i+1:n
            push!(vals, pfdmatrix[j,i])
        end
    end
    return vals
end

