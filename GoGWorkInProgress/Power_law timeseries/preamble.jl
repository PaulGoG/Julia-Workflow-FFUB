#Initialise program startup
#####
using Dates, DataFrames, CSV
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")

#Add path to input data folder
cd(@__DIR__)

δ_optimization = true
triDimPlots = false
const noDigits = 2

datafileType = ".csv"
datafileName = "10628"
datafileTime = "hours"
#datafileCol = 11
datafileCol = 9
#datafileQuantity = "atmospheric pressure (millibars)"
datafileQuantity = "wind speed (km/h)"

δ = 5.0