#Convert input data into desired format, standalone program

include("input.jl")
include("aux_func.jl")

using CSV, DataFrames, Tables

Edefdataname = "EXTRADEF.DSE"
dataEdef = CSV.read(Edefdataname, DataFrame; delim = ' ', ignorerepeated = true, skipto = 2, header = ["A", "Z", "Val_L", "Val_H"])

OutputEdef = zeros(length(dataEdef.A), 3)

for i in eachindex(dataEdef.A)
    A_H = Int(dataEdef.A[i])
    Z_H = Int(dataEdef.Z[i])
    Edef_L = dataEdef.Val_L[i]
    Edef_H = dataEdef.Val_H[i]
    A_L = A₀ - A_H
    Z_L = Z₀ - Z_H
    OutputEdef[i, 1] = A_H
    OutputEdef[i, 2] = Z_H
    OutputEdef[i, 3] = Edef_H
    if A_L != A_H 
        aux = [A_L Z_L Edef_L]
        OutputEdef = vcat(OutputEdef, aux)
    end
end

Output = Tables.table(OutputEdef)
CSV.write("Data.txt", Output, writeheader=false, newline = "\r\n")