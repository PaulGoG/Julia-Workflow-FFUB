if datafileType == ".CSV" || datafileType == ".csv"
    println("*reading $(string(datafileName, datafileType))...")
    cd("inputData/")

    #Data procured from the internet, *datafileCol* specifies which column has data of interest
    filedata = CSV.read(string(datafileName, datafileType), DataFrame; delim = ",")
    
    datapoints = DataFrame(Time = [i for i in first(axes(filedata))], Data = filedata[!, datafileCol])
    deleteat!(datapoints, findall(x -> ismissing(x), datapoints.Data))
    cd(@__DIR__)
elseif datafileType == ".DAT" || datafileType == ".dat"
    println("*reading $(string(datafileName, datafileType))...")
    cd("inputData/")

    #Manually curated data with no missing entries; data starts on second line and columns are delimited by a space
    datapoints = CSV.read(string(datafileName, datafileType), DataFrame;
    delim = " ", ignorerepeated = true, header = ["Time", "Data"], skipto = 2)
    cd(@__DIR__)
elseif datafileType == "dummy"
  #=First point of data is sampled from Gaussian of average 150 and variance 50
    The following point is sampled either from a Gaussian of average = lastDataPoint and variance = 3*δ or
    by another sampling of a Gaussian of average 150 and variance 50 (10% threshold)=#
    using Distributions
    println("*reading dummy Gaussian data...")
    datapoints = DataFrame(Time = Int[], Data = Float64[])
    push!(datapoints.Data, rand(Normal(150, 50)))
    push!(datapoints.Time, 1)
    for t in 2:Int(1e6)
        threshold = rand()
        if threshold <= 0.9
            push!(datapoints.Data, rand(Normal(datapoints.Data[t-1], 3*δ)))
        else
            push!(datapoints.Data, rand(Normal(150, 50)))
        end
        push!(datapoints.Time, t)
    end
else
    error("invalid input data type!")
end


datapoints.Time .= datapoints.Time .- datapoints.Time[1] .+ 1
datapoints.Data .= round.(datapoints.Data, digits = noDigits)