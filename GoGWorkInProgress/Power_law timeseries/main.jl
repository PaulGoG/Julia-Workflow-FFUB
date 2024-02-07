using Distributions, Plots, DataFrames, Optim

#Construction of test data for time series distribution

testData = rand(Normal(150, 50), Int(1e5))
testData = round.(testData, digits = 2)
testDays = 0:length(testData)
testδ = 5
testTimeseries = Float64[]

for i in eachindex(testData)
    searchableData = testData[i:length(testData)]
    index_delta = findfirst(x -> x >= first(searchableData) + testδ, searchableData)
    if !isnothing(index_delta)
        push!(testTimeseries, testDays[index_delta + i] - testDays[i])
    end
end

#Construct PDF and CDF of test timeseries

S = DataFrame(time = Float64[], PDF = Float64[], CDF = Float64[])
sort!(testTimeseries)

while(isassigned(testTimeseries, 1))
    indexes = findall(x -> x == testTimeseries[1], testTimeseries)
    push!(S.time, testTimeseries[1])
    push!(S.PDF, length(indexes))
    push!(S.CDF, sum(S.PDF))
    deleteat!(testTimeseries, indexes)
end

normalizationFactor = last(S.CDF)
S.PDF ./= normalizationFactor
S.CDF ./= normalizationFactor

#Construct power law model function and Kolmogorov-Smirnov loss function
function P_CDF(x, parameters)
    α, xₘ = parameters
    return (x/xₘ)^(1 - α)
end
function lossFunction(parameters, xdata, ydata) 
    xₘ = last(parameters)
    xmodel = xdata[xdata .>= xₘ]
    ymodel = [P_CDF(xmodel[i], parameters) for i in eachindex(xmodel)]
    ydataPlaceholder = ydata[xdata .>= xₘ]
    D = abs.(ymodel .- ydataPlaceholder)
    return maximum(D)
end

#Run Optimizer model 
xdata = S.time
ydata = S.CDF
initial_params = [2.0, minimum(xdata)]

optimizer_result = optimize(parameters -> lossFunction(parameters, xdata, ydata), initial_params)
α, xₘ = Optim.minimizer(optimizer_result)

P_PDF(x, α, xₘ) = (α - 1)/xₘ * (xₘ/x)^α

scatter(S.time, S.PDF, yscale = :log10)
plot!(S.time, P_PDF.(S.time, α, xₘ))

