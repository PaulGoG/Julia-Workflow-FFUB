#Blocks of code used only once for data conversion || Deprecated garbage code.

cd(@__DIR__)
cd("input_data/")

using CSV, DataFrames, Tables, Plots

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
#=
#Check s-wave neutron strength function from RIPL vs dse_eq_solvers
rawdatafile_name = "resonances0.dat"
rawdatafile_header = ["Z", "Symbol", "A", "Io", "Bn", "D0", "dD", "S0", "dS", "Gg", "dG", "Com"]
rawdatafile_firstline = 5
rawdatafile_delim = ' '
rawdatafile = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, skipto = rawdatafile_firstline, header = rawdatafile_header)
rawdatafile.S0 .*= 1e-4

RIPL_data = DataFrame(A = rawdatafile.A[rawdatafile.S0 .> 0], Z = rawdatafile.Z[rawdatafile.S0 .> 0], S_0 = rawdatafile.S0[rawdatafile.S0 .> 0])
A_data = Float64[]
S0_A_data = Float64[]
for A in sort(unique(RIPL_data.A))
    Denominator = 0.0
    Numerator = 0.0
    for Z in RIPL_data.Z[RIPL_data.A .== A]
        Numerator += RIPL_data.S_0[(RIPL_data.A .== A) .& ((RIPL_data.Z .== Z))][1]
        Denominator += 1.0
    end
    if Denominator > 0
        push!(S0_A_data, Numerator/Denominator)
        push!(A_data, A)
    end
end

scatter(A_data, S0_A_data)
S0_Calc = Strength_function_S₀.(A_data)
plot!(A_data, S0_Calc, yscale=:identity, size = (1000, 1000))
=#
#####

