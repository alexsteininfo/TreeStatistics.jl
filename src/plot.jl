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
    ylabel = cumulative ? "Cumulative number of mutations" : "Number of mutations"

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
    ylabel --> ylabel
    xlabel --> "VAF"
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
    ylabel = cumulative ? "Cumulative numberof mutations" : "Number of mutations"

    @series begin
        ylabel --> ylabel
        xlabel --> "Inverse VAF"
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
            ylabel --> ylabel
            xlabel --> "Inverse VAF"
            seriestype --> :line
            linecolor --> :red
            linestyle --> :dash
            x, y
        end
    end

    ylabel --> ylabel
    xlabel --> "Inverse VAF"
    fvals = [df[end, :VAF], 2/(1/df[1, :VAF]+ 1/df[end, :VAF]), df[1, :VAF]] 
    xticks --> (1 ./ fvals, map(x -> "1/" * x, string.(round.(fvals, digits = 2))))
    legend --> false
    grid --> false
    ()
end


@userplot PlotPopulation

@recipe function f(pp::PlotPopulation)
    for sresult in pp.args
        @series begin
            ylabel --> "Population size"
            xlabel --> "Time"
            seriestype --> :line
            sresult.output.tvec, sresult.output.Nvec
        end
    end
end