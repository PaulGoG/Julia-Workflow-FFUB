#Initialise program startup
#####
using Dates
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*loading input parameters and Julia libraries")

#Load Julia packages for data manipulation
using DataFrames, CSV

#Add path to input data folder
cd(@__DIR__)

#Check for input folder existance
if !isdir("input_data/")
    mkdir("input_data/")
end
cd("input_data/")

#Define main struct objects
abstract type AbstractDistribution end
struct Distribution{T1 <: Vector{Int}, T2 <: Vector{Float64}} <: AbstractDistribution
    A::T1
    Z::T1
    TKE::T2
    No_Sequence::T1
    Value
    σ::T2
end
struct Distribution_unidym{T <: Vector{Float64}} <: AbstractDistribution
    Argument
    Value::T
    σ::T
end