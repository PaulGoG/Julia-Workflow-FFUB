#Construction power law model function and Kolmogorov-Smirnov loss function
using Optim

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
    return maximum(d_KS)
end