"""
In this code we are going to simulate atmospheric tritium emission and dispersion
for normal nuclear reactor working conditions.
"""
# Main code body

using Plots
using Trapz
using DataFrames
using CSV

include("RawDataRead.jl")


include("TestDependency.jl")
a="A"
Selector_Constanta(a)




