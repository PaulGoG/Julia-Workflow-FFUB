#Initialise program startup
#####
using Dates
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")

#Add path to input data folder
cd(@__DIR__)

δ_optimization = false
datafileTimeHourly = true
triDimPlots = false

datafileType = ".csv"
datafileName = "10628"
datafileCols = 9
noDecimals = 3