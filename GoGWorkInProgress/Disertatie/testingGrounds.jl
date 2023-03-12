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

RIPL3_data = DataFrame(A = rawdatafile.A[rawdatafile.S0 .> 0], Z = rawdatafile.Z[rawdatafile.S0 .> 0], S_0 = rawdatafile.S0[rawdatafile.S0 .> 0])

function Strength_function_new(A)
    if A <= 23
        return 2.5e-5
    elseif A > 23 && A <= 55
        return 1.296875e-5 *A - 0.00027328125
    elseif A > 55 && A <= 69
        return -2.142857142857143e-5 *A + 0.0016185714285714289
    elseif A > 69 && A <= 79
        return 0.00014
    elseif A > 79 && A <= 90
        return -7.818181818181818e-6 *A + 0.0007576363636363635
    elseif A > 90 && A <= 110
        return 5.4e-5
    elseif A > 110 && A <= 120
        return -4.4e-6 *A + 0.000538 
    elseif A > 120 && A <= 135
        return 6.666666666666667e-6 *A - 0.00079
    elseif A > 135 && A <= 146
        return 2.4545454545454552e-5 *A - 0.0032036363636363647
    elseif A > 146 && A <= 152
        return 0.00038
    elseif A > 152 && A <= 158
        return -3.833333333333334e-5 *A + 0.006206666666666668
    elseif A > 158 && A <= 173
        return 0.00015
    elseif A > 173 && A <= 186
        return 1.0769230769230771e-5 *A - 0.0017130769230769235
    elseif A > 186 && A <= 208
        return -8.636363636363637e-6 *A + 0.0018963636363636366
    elseif A > 208
        return 1e-4
    end
end

function Strength_function_old(A)
    if A <= 70
        return 7e-5
    elseif A > 70 && A <= 86
        return 1e-4 *(A*1.875e-2 - 6.125e-1)
    elseif A > 86 && A <= 111
        return 1e-4 
    elseif A > 111 && A <= 121
        return 1e-4 *(-A*2.857e-2 + 4.1714)
    elseif A > 121 && A <= 140
        return 1e-4 
    elseif A > 140 && A <= 144
        return 1e-4 *(A*7.5e-2 - 9.8)
    elseif A > 144
        return 1e-4
    end
end

function slope_intercept(x_1, y_1, x_2, y_2)
    a = (y_1-y_2)/(x_1-x_2)
    b = y_1 - a*x_1
    return a, b
end

slope_intercept(135, 0.00010999999999999996, 146, 380e-6)

S0_old = Strength_function_old.(sort(unique(RIPL3_data.A)))
S0_new = Strength_function_new.(sort(unique(RIPL3_data.A)))

pltlog = scatter(RIPL3_data.A, RIPL3_data.S_0, yscale=:log10, size = (1280, 720), color = :black, xlabel = "A", ylabel = "S₀", framestyle = :box, minorgrid = true, title = "logarithmic scale")
pltlog = scatter!(pltlog, sort(unique(RIPL3_data.A)), S0_old, color = :red)
pltlog = plot!(pltlog, sort(unique(RIPL3_data.A)), S0_new, color = :blue)

pltlin = scatter(RIPL3_data.A, RIPL3_data.S_0, yscale=:identity, size = (1280, 720), color = :black, xlabel = "A", ylabel = "S₀", framestyle = :box, minorgrid = true, title = "linear scale")
pltlin = scatter!(pltlin, sort(unique(RIPL3_data.A)), S0_old, color = :red)
pltlin = plot!(pltlin, sort(unique(RIPL3_data.A)), S0_new, color = :blue)

=#
#####

