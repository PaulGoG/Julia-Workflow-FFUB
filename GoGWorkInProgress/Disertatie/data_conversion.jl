#Convert input data into desired format, standalone program, program abstraction kept to a minimum!

include("input.jl")
include("aux_func.jl")

using Tables

#=
!DATA-SPECIFIC CODE!
Blocks of code specific to one type of file conversion will be commented away
=#
#####
#Convert L&H stacked values provided on H domain to single-file values on the whole A domain
#=
rawdatafile_name = "EXTRADEF.DSE"
rawdatafile_header = ["A", "Z", "Val_L", "Val_H"]
rawdatafile_firstline = 2
rawdatafile_delim = ' '
rawdatafile = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, skipto = rawdatafile_firstline, header = rawdatafile_header)

Output = zeros(nrow(rawdatafile), 3)
for i in eachindex(Output[:, 1])
    A_H = rawdatafile[!, 1][i]
    Z_H = rawdatafile[!, 2][i]
    Value_L = rawdatafile[!, 3][i]
    Value_H = rawdatafile[!, 4][i]
    A_L = A₀ - A_H
    Z_L = Z₀ - Z_H
    Output[i, 1] = A_H
    Output[i, 2] = Z_H
    Output[i, 3] = Value_H
    if A_L != A_H 
        Output = vcat(Output, [A_L Z_L Value_L])
    end
end

Output = Tables.table(Output;  header=[:A, :Z, :Value])
newdatafile_name = "Converted_$rawdatafile_name"
CSV.write(newdatafile_name, Output, delim="   ")
=#

#Concatenate data for ΔZ(A) & rms(A) in a single file
#=
rawdatafile_name = "U5DELTAZ.WAH"
rawdatafile_header = ["A", "Val"]
rawdatafile_firstline = 2
rawdatafile_delim = ' '
rawdatafile_1 = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, skipto = rawdatafile_firstline, header = rawdatafile_header)

rawdatafile_name = "U5RMS.WAH"
rawdatafile_2 = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, skipto = rawdatafile_firstline, header = rawdatafile_header)

if nrow(rawdatafile_1) == nrow(rawdatafile_2)
    Output =  hcat(rawdatafile_1[!, 1], rawdatafile_1[!, 2], rawdatafile_2[!, 2])
    Output = Tables.table(Output;  header=[:A, :ΔZ, :rms])
    newdatafile_name = "Converted_deltaZA_rmsA.CSV"
    CSV.write(newdatafile_name, Output, delim="   ")
else error("raw data files row sizes do not match!")
end
=#
#####