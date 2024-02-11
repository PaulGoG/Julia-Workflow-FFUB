#Nested optimization of best δ for the input dataset that minimizez the Kolmogorov-Smirnov statistic
function lossDelta(datapoints::DataFrame, δ)
    #If δ number of digits is unconstrained it will tend towards 0!
    if !iszero(floor(only(δ), digits = noDigits))
        println("δ = $(only(δ))")
        KS = powerlawOptimize(wtsProbabilities(waiting_time_series(datapoints, only(δ))))[3]
        if !isnan(KS)
            return KS
        else
            lossDelta(datapoints, δ/10)*100
        end
    else 
        println("δ value too close to 0, retrying search...")
        return Inf
    end
end
function optimumDelta(datapoints::DataFrame, starting_δ::Float64)
    println("*searching for optimum δ value...\n")
    initial_parameter = [starting_δ]
    result = optimize(parameter -> lossDelta(datapoints, round.(abs.(parameter), digits = noDigits)), initial_parameter)
    δ_optim = round(abs(only(Optim.minimizer(result))), digits = noDigits)
    println("\n*!Warning!: If δ value is too many orders of magnitute lower than the average input data provided,\n     it may indicate non-FreeScale statistical behaviour!\n")
    println("*Optimum δ value found, printing status log:\n")
    display(result)
    return δ_optim
end