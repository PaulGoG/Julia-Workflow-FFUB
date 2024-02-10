#Construction of waiting time series between Data and Data + δ at each point
#Brute force O(n^2) solution
function waiting_time_series(datapoints::DataFrame, δ)
    println("*computing time series...")
    WTS = Int[]
    for index in eachindex(datapoints.Data)
        searchableData = @view datapoints.Data[(index + 1):end]
        index_delta = findfirst(x -> x >= datapoints.Data[index] + δ, searchableData)
        if !isnothing(index_delta)
            push!(WTS, datapoints.Time[index + index_delta] - datapoints.Time[index])
        end
    end
    return sort!(WTS)
end
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
function wtsProbabilities(WTS::Array{Int})
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
    return timeseries
end