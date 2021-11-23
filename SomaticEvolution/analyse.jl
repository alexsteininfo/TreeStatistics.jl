function fitinverse(df::DataFrame, fmin, fmax)
    df = df[fmin .<= df.VAF .< fmax, :]
    #remove clonal mutations from the cumultive frequency
    df[!, :cumfreq] .= df[!, :cumfreq] .- df[end, :cumfreq] .+ df[end, :freq]
    #fit the data to y = m(x - 1/fmax), enforcing the intercept y(1/fmax) = 0
    #y = M(f), x = 1/f and the gradient m = μ/β
    df.invVAF = 1 ./df[!, :VAF] .- 1/fmax
    lmfit = fit(LinearModel, @formula(cumfreq ~ invVAF + 0), df)
    df.prediction = predict(lmfit)
    df.residuals = residuals(lmfit)
    df[!,:invVAF] .+= 1/fmax
    return df, coef(lmfit)[1], r2(lmfit)  
end

function fitinverse(VAF, fmin, fmax; xstep=0.001)
    VAF = filter(x-> fmin <= x <= fmax, VAF)
    df = gethist(VAF, xmin = fmin, xmax = fmax, xstep = xstep)
    #fit the data to y = m(x - 1/fmax), enforcing the intercept y(1/fmax) = 0
    #y = M(f), x = 1/f and the gradient m = μ/β
    df.invVAF = 1 ./df[!, :VAF] .- 1/fmax
    lmfit = fit(LinearModel, @formula(cumfreq ~ invVAF + 0), df)
    df.prediction = predict(lmfit)
    df.residuals = residuals(lmfit)
    df[!,:invVAF] .+= 1/fmax
    return df, coef(lmfit)[1], r2(lmfit)  
end