# Simulam dispersia atmosferica a gazului de tritiu care iese la
# cosul unui reactor nuclear de tip CANDU

using Plots
using Trapz
using DataFrames
using CSV

include("Constante.jl")
include("CitireDate.jl")