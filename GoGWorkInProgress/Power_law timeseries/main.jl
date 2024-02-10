include("preamble.jl")
include("readData.jl")
include("buildTimeseries.jl")
include("parameterOptimization.jl")

if δ_optimization
    include("optimizeDelta.jl")
else
    println("Input value for δ = ")
    δ = parse(Float64, readline())

    timeseries = wtsProbabilities(waiting_time_series(datapoints, δ))

    xData = timeseries.time
    yData = timeseries.cdf
    initial_params = [2.0, minimum(xData)]

    println("*parameter optimization inbound...")
    optimizer_result = optimize(parameters -> lossKolmogorovSmirnov(parameters, xData, yData), initial_params)
    display(optimizer_result)

    parameters = Optim.minimizer(optimizer_result)
    d_KS = Optim.minimum(optimizer_result)
    α, xₘ = parameters

    P_PDF(x, α, xₘ) = (α - 1)/xₘ * (xₘ/x)^α
    xFScale = @view xData[xData .>= xₘ]
end

include("plotsfigs.jl")

println("*end program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*program execution succesful!")