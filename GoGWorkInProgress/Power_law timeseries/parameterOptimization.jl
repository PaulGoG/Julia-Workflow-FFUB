#Construction power law model function and Kolmogorov-Smirnov loss function
using Optim
#α is the first parameter and xₘ is the second parameter
function Freescale_CDF(x, parameters)
    α, xₘ = parameters
    return (x/xₘ)^(1 - α)
end
function Pareto_Tsallis_CDF(x, parameters)
    α, λ = parameters
    return 1 - (1 + x/λ)^(-α)
end
function lossFS(parameters, xData, yData)
    α = first(parameters)
    xₘ = last(parameters)
    if α > 0 && xₘ > 0
        xFScale = @view xData[xData .>= xₘ]
        yFScale = [Freescale_CDF(xFScale[i], parameters) for i in eachindex(xFScale)]
        yDataPlaceholder = @view yData[xData .>= xₘ]
        d_KS = abs.(yFScale .- yDataPlaceholder)
        if !isempty(d_KS)
            return maximum(d_KS)
        else 
            return Inf
        end
    else
        return Inf
    end
end
function lossPT(parameters, xData, yData)
    α = first(parameters)
    λ = last(parameters)
    if α > 0 && λ > 0
        yPT = [Pareto_Tsallis_CDF(xData[i], parameters) for i in eachindex(xData)]
        d_KS = abs.(yPT .- yData)
        if !isempty(d_KS)
            return maximum(d_KS)
        else 
            return Inf
        end
    else
        return Inf
    end
end
function lossKolmogorovSmirnov(parameters, xData, yData, distribution_type)
    if distribution_type == "FS"
        lossFS(parameters, xData, yData)
    elseif distribution_type == "PT"
        lossPT(parameters, xData, yData)
    else 
        error("Invalid distribution_type = $(distribution_type)")
    end
end
function powerlawOptimize(timeseries::DataFrame, distribution_type)
    if !isempty(timeseries)
        if distribution_type == "FS"
            initial_params = [2.0, minimum(timeseries.time)]
            optimizer_result = optimize(parameters -> lossKolmogorovSmirnov(parameters, timeseries.time, timeseries.cdf, distribution_type), initial_params)
            α, xₘ = Optim.minimizer(optimizer_result)
            d_KS = Optim.minimum(optimizer_result)
            return α, xₘ, d_KS
        elseif distribution_type == "PT"
            initial_params = [2.0, sum(timeseries.time)/length(timeseries.time)]
            optimizer_result = optimize(parameters -> lossKolmogorovSmirnov(parameters, timeseries.time, timeseries.cdf, distribution_type), initial_params)
            α, λ = Optim.minimizer(optimizer_result)
            d_KS = Optim.minimum(optimizer_result)
            return α, λ, d_KS
        else 
            error("Invalid distribution_type = $(distribution_type)")
        end
    else
        return NaN, NaN, NaN
    end
end