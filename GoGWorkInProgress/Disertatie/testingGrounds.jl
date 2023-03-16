#Blocks of code used only once for data conversion || Deprecated garbage code.

cd(@__DIR__)
cd("input_data/")
using CSV, DataFrames, Plots

#Convert L&H stacked values provided on H domain to single-file values on the whole A domain
#=
rawdatafile_name = "CF2YATKE.VES"
rawdatafile_header = ["A", "TKE", "Val"]
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

aux_col = [NaN for i in 1:nrow(rawdatafile)]

Output = DataFrame(A = rawdatafile.A, TKE = rawdatafile.TKE, Value = rawdatafile.Val, σ = aux_col)
newdatafile_name = "CF2YATKE.VES"
CSV.write(newdatafile_name, Output, delim=' ')
=#

#Concatenate data for ΔZ(A) & rms(A) in a single file
#=
rawdatafile_name = "CFDELTAZ.WAH"
rawdatafile_header = ["A", "Val"]
rawdatafile_firstline = 2
rawdatafile_delim = ' '
rawdatafile_1 = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, skipto = rawdatafile_firstline, header = rawdatafile_header)

rawdatafile_name = "CFRMS.WAH"
rawdatafile_2 = CSV.read(rawdatafile_name, DataFrame; delim = rawdatafile_delim, ignorerepeated = true, skipto = rawdatafile_firstline, header = rawdatafile_header)

Output =  DataFrame(A = rawdatafile_1.A, ΔZ = rawdatafile_1.Val, rms = rawdatafile_2.Val)
newdatafile_name = "DeltaZ_rms_A.CF2"
CSV.write(newdatafile_name, Output, delim=' ')

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
slope_intercept(150, 0.9, 167, 0.68)

S0_old = Strength_function_old.(sort(unique(RIPL3_data.A)))
S0_new = Strength_function_new.(sort(unique(RIPL3_data.A)))

pltlog = scatter(
    RIPL3_data.A, RIPL3_data.S_0, 
    yscale=:log10, size = (1000, 900), color = :black, xlabel = "A", ylabel = "S₀", 
    framestyle = :box, minorgrid = true, title = "Log10 scale", 
    xticks = (20:10:300), label="RIPL3 data"
    )
pltlog = scatter!(pltlog, sort(unique(RIPL3_data.A)), S0_old, color = :red, label="S₀ old")
pltlog = plot!(pltlog, sort(unique(RIPL3_data.A)), S0_new, color = :blue, label="S₀ new")
vline!(pltlog, [76, 160], color = :green, label="")

pltlin = scatter(
    RIPL3_data.A, RIPL3_data.S_0, 
    yscale=:identity, size = (1000, 900), color = :black, xlabel = "A", ylabel = "S₀", 
    framestyle = :box, minorgrid = true, title = "Lin scale",
    xticks = (20:10:300), label="RIPL3 data"
    )
pltlin = scatter!(pltlin, sort(unique(RIPL3_data.A)), S0_old, color = :red, label="S₀ old")
pltlin = plot!(pltlin, sort(unique(RIPL3_data.A)), S0_new, color = :blue, label="S₀ new")
vline!(pltlin, [76, 160], color = :green, label="")
=#
#Generate ΔE_def for given β₀ & β_scission parametrization
#=
dβ₀ = CSV.read("B2MOLLER.ANA", DataFrame; delim=' ', ignorerepeated=true, header=["Z", "A", "Value"], skipto = 2)
function Linear_function(x, a, b)
    return a*x + b
end
function Beta_sciss(Z)
    if Z < 28
        return NaN
    elseif Z >= 28 && Z <= 41
        a = 0.58/13
        b = -28*a
        return Linear_function(Z, a, b)
    elseif Z >= 42 && Z <= 44
        return 0.58
    elseif Z >= 45 && Z <= 50
        a = -0.58/6
        b = -50*a
        return Linear_function(Z, a, b)
    elseif Z >= 51 && Z <= 70
        a = 0.6/15
        b = -50*a
        return Linear_function(Z, a, b)
    else
        return NaN
    end
end
function E_LDM(β, A, Z)
    η = (A - 2*Z)/A
    χ = 1 - 1.7826 * η^2
    α² = (5 * β^2)/(4*π)
    X_νs = -χ*(15.4941*A - 17.9439*A^(2/3) * (1 + 0.4*α²))
    X_e = Z^2 *(0.7053 * (1 - 0.2*α²)/(A^(1/3)) - 1.1529/A)
    return X_νs + X_e
end
function ΔE_def(dβ₀)
    ΔE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in sort(unique(dβ₀.A))
        for Z in sort(unique(dβ₀.Z[dβ₀.A .== A]))
            β_sciss = Beta_sciss(Z)
            β_0 = dβ₀.Value[(dβ₀.A .== A) .& (dβ₀.Z .== Z)][1]
            if !isnan(β_sciss)
                E_LDM_0 = E_LDM(0, A, Z)
                E_def_sciz = E_LDM(β_sciss, A, Z) - E_LDM_0
                E_def_at = E_LDM(β_0, A, Z) - E_LDM_0
                δE_def = abs(E_def_sciz - E_def_at)
                push!(ΔE.A, A)
                push!(ΔE.Z, Z)
                push!(ΔE.Value, δE_def)
            end
        end
    end
    return ΔE
end
ΔE_deformation = ΔE_def(dβ₀)
Output =  DataFrame(A = ΔE_deformation.A, Z = ΔE_deformation.Z, ΔE_def = ΔE_deformation.Value)
newdatafile_name = "EXTRADEF.U5"
CSV.write(newdatafile_name, Output, delim=' ')
=#
#####

