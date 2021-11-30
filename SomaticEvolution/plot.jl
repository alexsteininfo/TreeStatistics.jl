@userplot struct PlotVAF{T<:Tuple{Simulation}}
    args::T
end

@recipe function f(pv::PlotVAF; sampled = true, cumulative = false, xstep = 0.01)
    if sampled    
        df = gethist(pv.args[1].sampled.VAF, fstep = xstep)
    else
        df = gethist(pv.args[1].output.trueVAF, fstep = xstep)
    end
    VAF = df[:,:VAF]
    freq = cumulative ? df[:,:cumfreq] : df[:,:freq]
    ylabel = cumulative ? "Cumulative number of mutations" : "Number of mutations"

    @series begin
        seriestype --> :bar
        linecolor --> :white
        linewidth --> 0.2
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


@recipe function f(piv::PlotInverseVAF; sampled = true, fmin = 0.12, fmax = 0.24, 
                    fstep = 0.001, fitcoef = nothing, cumulative = true, addtext = true, 
                    rsq = nothing, plottrue = false)
    if sampled 
        df = gethist(piv.args[1].sampled.VAF, fmin = fmin, fmax = fmax, fstep = fstep) 
    else 
        df = gethist(piv.args[1].output.trueVAF, fmin = fmin, fmax = fmax, fstep = fstep)
    end
    VAF = df[!,:VAF]
    freq = cumulative ? df[!,:cumfreq] : df[!,:freq]
    ylabel = cumulative ? "Cumulative number of mutations" : "Number of mutations"

    # if addtext
    #     β = 1 - piv.args[1].input.siminput.d/piv.args[1].input.siminput.b
    #     μ = piv.args[1].input.siminput.μ
    #     text = "μ/β_pred = $fitcoef"
    #     text *= "\nμ/β_true = $(μ/β)"
    #     text *= rsq !== nothing ? "\nr^2 = $rsq" : ""
    #     annotation = (1/0.16, 500, text)
    # else annotation = [] end

    @series begin
        ylabel --> cumulative ? "Cumulative number of mutations" : "Number of mutations"
        xlabel --> "Inverse VAF"
        seriestype --> :line
        # annotationhalign --> :right
        # annotationvalign --> :bottom
        # annotation := annotation
        1 ./ VAF, freq
    end

    if fitcoef !== nothing
        @series begin
            ylabel --> cumulative ? "Cumulative number of mutations" : "Number of mutations"
            xlabel --> "Inverse VAF"
            seriestype --> :line
            linecolor --> :red
            linestyle --> :dash
            1 ./ VAF, fitcoef .* (1 ./ VAF .- 1/fmax)
        end
    end

    # if plottrue
    #     β = 1 - piv.args[1].input.siminput.d/piv.args[1].input.siminput.b
    #     μ = piv.args[1].input.siminput.μ
    #     @series begin
    #         ylabel --> cumulative ? "Cumulative number of mutations" : "Number of mutations"
    #         xlabel --> "Inverse VAF"
    #         seriestype --> :line
    #         linecolor --> :green
    #         linestyle --> :dash
    #         1 ./ VAF, μ/β .* (1 ./ VAF .- 1/fmax)
    #     end
    # end

    ylabel --> ylabel
    xlabel --> "Inverse VAF"
    fvals = [fmax, 2/(1/fmin + 1/fmax), fmin] 
    xticks --> (1 ./ fvals, map(x -> "1/" * x, string.(round.(fvals, digits = 2))))
    legend --> false
    grid --> false
    ()
end
