using DataFrames, CSV
cd("inputData/")
println("*reading $(string(datafileName, datafileType))...")

if datafileType == ".CSV" || datafileType == ".csv"
    filedata = CSV.read(string(datafileName, datafileType), DataFrame; delim = ",")
    
    datapoints = DataFrame(Time = [i for i in first(axes(filedata))], Data = filedata[!, datafileCols])
    deleteat!(datapoints, findall(x -> ismissing(x), datapoints.Data))
elseif datafileType == ".DAT" || datafileType == ".dat"
    datapoints = CSV.read(string(datafileName, datafileType), DataFrame;
    delim = " ", ignorerepeated = true, header = ["Time", "Data"], skipto = 2)
else
    error("invalid input data type!")
end

cd(@__DIR__)

if !datafileTimeHourly
    println("Time multiplier for hourly = ")
    timeMultiplier = parse(Int, readline())
    datapoints.Time .*= timeMultiplier
end

datapoints.Time .-= datapoints.Time[1] + 1
datapoints.Data .= round.(datapoints.Data, digits = noDecimals)