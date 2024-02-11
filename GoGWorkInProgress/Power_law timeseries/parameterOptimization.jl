#Construction power law model function and Kolmogorov-Smirnov loss function
using Optim
#α is the first parameter and xₘ is the second parameter
function P_CDF(x, parameters)
    α, xₘ = parameters
    return (x/xₘ)^(1 - α)
end
function lossKolmogorovSmirnov(parameters, xData, yData)
    xₘ = last(parameters)
    xFScale = @view xData[xData .>= xₘ]
    yFScale = [P_CDF(xFScale[i], parameters) for i in eachindex(xFScale)]
    yDataPlaceholder = @view yData[xData .>= xₘ]
    d_KS = abs.(yFScale .- yDataPlaceholder)
    if !isempty(d_KS)
        return maximum(d_KS)
    else 
        return Inf
    end
end
function powerlawOptimize(timeseries::DataFrame)
    if !isempty(timeseries)
        initial_params = [2.0, minimum(timeseries.time)]
        optimizer_result = optimize(parameters -> lossKolmogorovSmirnov(parameters, timeseries.time, timeseries.cdf), initial_params)
        α, xₘ = Optim.minimizer(optimizer_result)
        d_KS = Optim.minimum(optimizer_result)
        return α, xₘ, d_KS
    else
        return NaN, NaN, NaN
    end
end