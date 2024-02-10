using Distributions, DataFrames, Optim, Plots

#Construction of dummy input data for the optimiser
#First point of data is sampled from Gaussian of average 150 and variance 50
#The following point is sampled either from a Gaussian of average = lastDataPoint and variance = 3*δ or
#by another sampling of a Gaussian of average 150 and variance 50 (5% threshold)

datapointsLength = Int(1e5)
δ = 10.0
noDecimals = 2

datapoints = DataFrame(Data = Float64[], Time = Int[])
push!(datapoints.Data, rand(Normal(150, 50)))
push!(datapoints.Time, 1)
for t in 2:datapointsLength
    threshold = rand()
    if threshold <= 0.95
        push!(datapoints.Data, rand(Normal(datapoints.Data[t-1], 3*δ)))
    else
        push!(datapoints.Data, rand(Normal(150, 50)))
    end
    push!(datapoints.Time, t)
end
datapoints.Data .= round.(datapoints.Data, digits = noDecimals)

#Construction of waiting time series between Data and Data + δ at each point
#Brute force O(n^2) solution
WTS = Int[]
for index in eachindex(datapoints.Data)
    searchableData = @view datapoints.Data[(index + 1):end]
    index_delta = findfirst(x -> x >= datapoints.Data[index] + δ, searchableData)
    if !isnothing(index_delta)
        push!(WTS, datapoints.Time[index + index_delta] - datapoints.Time[index])
    end
end
sort!(WTS)

#=
#Refined O(n) solution
#Work still in progress!
function waiting_time_series(datapoints::DataFrame, δ)
    indexStack = Int[]
    timeseries = Int[]
    for index in eachindex(datapoints.Data)
        while !isempty(indexStack) && datapoints.Data[index] >= datapoints.Data[last(indexStack)] + δ
            prev_index = pop!(indexStack)
            push!(timeseries, datapoints.Time[index] - datapoints.Time[prev_index])
        end
        push!(indexStack, index)
    end
    return sort!(timeseries)
end
=#

#Construction of PDF and CDF for the time series data
timeseries = DataFrame(time = Int[], pdf = Float64[], cdf = Float64[])

while(!isempty(WTS))
    indices = findall(x -> x == WTS[1], WTS)
    push!(timeseries.time, WTS[1])
    push!(timeseries.pdf, length(indices))
    push!(timeseries.cdf, sum(timeseries.pdf))
    deleteat!(WTS, indices)
end
normalizationFactor = last(timeseries.cdf)
timeseries.pdf ./= normalizationFactor
timeseries.cdf ./= normalizationFactor

#deleteat!(timeseries, findall(x -> x <= 1e-4, timeseries.pdf))

#Construction power law model function and Kolmogorov-Smirnov loss function
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

#Run optimizer model on data and basic plot
xData = timeseries.time
yData = timeseries.cdf
initial_params = [2.0, minimum(xData)]

optimizer_result = optimize(parameters -> lossKolmogorovSmirnov(parameters, xData, yData), initial_params)
parameters = Optim.minimizer(optimizer_result)
d_KS = Optim.minimum(optimizer_result)
α, xₘ = parameters

P_PDF(x, α, xₘ) = (α - 1)/xₘ * (xₘ/x)^α
xFScale = @view xData[xData .>= xₘ]

scatter(xData, timeseries.pdf, yscale = :log10)
plot!(xFScale, P_PDF.(xFScale, α, xₘ))
yticks!([10.0^i for i in -10:10])