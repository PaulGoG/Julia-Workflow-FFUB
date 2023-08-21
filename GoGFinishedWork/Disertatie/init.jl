#Initialise program startup
#####
using Dates
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*loading input parameters and Julia libraries")

#Add path to input data folder
cd(@__DIR__)

#Check for input folder existance
if !isdir("input_data/")
    mkdir("input_data/")
end
cd("input_data/")