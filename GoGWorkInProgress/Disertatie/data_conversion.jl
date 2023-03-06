#Convert input data into desired format, standalone program, program abstraction kept to a minimum!

cd(@__DIR__)
cd("input_data/")

using CSV, DataFrames, Tables

#=
!DATA-SPECIFIC CODE!
Blocks of code specific to one type of file conversion
=#
#####
#Convert L&H stacked values provided on H domain to single-file values on the whole A domain
#=
rawdatafile_name = "EXTRADF5.DSE"
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
newdatafile_name = "EXTRADEF.IN"
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
    newdatafile_name = "DeltaZA_rmsA.U5"
    CSV.write(newdatafile_name, Output, delim="   ")
else error("raw data files row sizes do not match!")
end
=#
#####
#=
    #Sorts data in ascending order for mass number
    length_Y = length(y_A_Z_TKE.Value)
    aux_A = sort(y_A_Z_TKE.A)
    aux_Z = [y_A_Z_TKE.Z[y_A_Z_TKE.A .== A] for A in unique(aux_A)]
    aux_Z = reduce(vcat, aux_Z)
    aux_TKE = zeros(length_Y)
    aux_Value = zeros(length_Y)
    aux_σ = zeros(length_Y)
    aux_index = 0
    for index_fragmdomain in eachindex(fragmdomain.A)
        for TKE in unique(y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== fragmdomain.A[index_fragmdomain]) .& (y_A_Z_TKE.Z .== fragmdomain.Z[index_fragmdomain])])
            aux_index += 1
            aux_TKE[aux_index] = y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== fragmdomain.A[index_fragmdomain]) .& (y_A_Z_TKE.Z .== fragmdomain.Z[index_fragmdomain]) .& (y_A_Z_TKE.TKE .== TKE)][1]
            aux_Value[aux_index] = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== fragmdomain.A[index_fragmdomain]) .& (y_A_Z_TKE.Z .== fragmdomain.Z[index_fragmdomain]) .& (y_A_Z_TKE.TKE .== TKE)][1]
            aux_σ[aux_index] = y_A_Z_TKE.σ[(y_A_Z_TKE.A .== fragmdomain.A[index_fragmdomain]) .& (y_A_Z_TKE.Z .== fragmdomain.Z[index_fragmdomain]) .& (y_A_Z_TKE.TKE .== TKE)][1]
        end
    end
    for index in 1:aux_index
        y_A_Z_TKE.A[index] = aux_A[index]
        y_A_Z_TKE.Z[index] = aux_Z[index]
        y_A_Z_TKE.TKE[index] = aux_TKE[index]
        y_A_Z_TKE.Value[index] = aux_Value[index]
        y_A_Z_TKE.σ[index] = aux_σ[index]
    end
=#