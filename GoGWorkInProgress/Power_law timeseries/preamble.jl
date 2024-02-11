#Initialise program startup
#####
using Dates, DataFrames, CSV
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")

#Add path to input data folder
cd(@__DIR__)

δ_optimization = true
triDimPlots = false
const noDigits = 3

datafileType = ".csv"
datafileName = "meteoBucurestiImh"
datafileTime = "days"
datafileCol = 10
datafileQuantity = "atmospheric pressure (millibars)"

δ = 1000.0