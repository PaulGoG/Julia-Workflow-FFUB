#Initialise program startup
#####
using Dates, DataFrames, CSV
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")

#Add path to input data folder
cd(@__DIR__)

δ_optimization = false
triDimPlots = false
const noDigits = 2

distribution_type = "PT"

datafileType = ".csv"
datafileName = "10628"
datafileTime = "hours"
datafileCol = 11
#datafileCol = 9
datafileQuantity = "atmospheric pressure (millibars)"
#datafileQuantity = "wind speed (m/s)"

δ = 50.0