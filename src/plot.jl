@userplot struct PlotVAF{T<:Tuple{Simulation}}
    args::T
end

@recipe function f(pv::PlotVAF; sampled = true, cumulative = false, fstep = 0.01)
    if sampled    
        df = gethist(pv.args[1].sampled.VAF, fstep = fstep)
    else
        df = gethist(pv.args[1].output.trueVAF, fstep = fstep)
    end
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

    if length(pv.args[1].output.subclones) > 0
        subclonefreq = [subclone.freq for subclone in pv.args[1].output.subclones]
        xint = subclonefreq ./2 * pv.args[1].input.cellularity

        @series begin
            seriestype --> :vline
            legend --> false
            fillcolor --> :darkred
            linewidth --> 2
            legend --> false
            grid --> false
            xint
        end
    end
    yguide --> yguide
    xguide --> "VAF"
    legend --> false
    grid --> false
    ()
end

@userplot struct PlotInverseVAF{T<:Tuple{Simulation}}
    args::T
end


@recipe function f(piv::PlotInverseVAF; sampled = true, fmin = 0.12, fmax = 0.24, fit_fmax=fmax,
                    fstep = 0.001, fitcoef = nothing, cumulative = true, dataseries = :line)
    if sampled 
        df = gethist(piv.args[1].sampled.VAF, fmin = fmin, fmax = fmax, fstep = fstep) 
    else 
        df = gethist(piv.args[1].output.trueVAF, fmin = fmin, fmax = fmax, fstep = fstep)
    end
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



@recipe function f(output::SimulationResult)
    @series begin
        yguide --> "Population size"
        xguide --> "Time"
        seriestype --> :line
        output.tvec, output.Nvec
    end
end

@recipe function f(outputvec::Vector{SimulationResult}; plottype=:popsize)
    if plottype == :modulesize
        yguide --> "Module size"
        xguide --> "Time"
        for output in outputvec
            @series begin
                seriestype --> :line
                legend --> false
                output.tvec, output.Nvec
            end
        end
    elseif plottype == :popsize
        yguide --> "Number of modules"
        xguide --> "Time"
        @series begin
            seriestype --> :line
            legend --> false
            newmoduletimes(outputvec), 1:length(outputvec)
        end
    # elseif plottype == :popsizebycell
    #     yguide --> "Number of cells"
    #     xguide --> "Time"
    #     @series begin
            
    #     end

    end
end

@userplot PairwisePlot

@recipe function f(pp::PairwisePlot)
    if length(pp.args) == 1
        z = pp.args[1]
        x, y = 1:size(z)[1], 1:size(z)[2]
    elseif length(pp.args) == 3
        x, y, z = pp.args
    end
    yguide --> "Module #"
    xguide --> "Module #"
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
