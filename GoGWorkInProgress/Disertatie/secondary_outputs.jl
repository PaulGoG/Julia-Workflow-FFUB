#=
This part of the program processes the main output so that
it can be compared with experimental data (averaging main output data over experimental yield distributions)
=#
#####
#Compute Y(A,Z,TKE) form experimental Y(A,TKE) data
function Process_yield_data(A_0, fragmdomain::Distribution, dY::DataFrame)
    y_A_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    #Completes fragmentation domain in case data provided only for HF
    if !isassigned(dY.A[dY.A .< A_0/2], 1)
        for A in sort(unique(dY.A), rev=true)
            for Z in fragmdomain.Z[fragmdomain.A .== A_0 - A]
                P_A_Z = fragmdomain.Value[(fragmdomain.A .== A_0 - A) .& (fragmdomain.Z .== Z)][1]
                for TKE in sort(unique(dY.TKE[dY.A .== A]))
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    if A != A_0 - A
                        push!(y_A_Z_TKE.A, A_0 - A)
                        push!(y_A_Z_TKE.Z, Z)
                        push!(y_A_Z_TKE.TKE, TKE)
                        push!(y_A_Z_TKE.Value, val)
                        push!(y_A_Z_TKE.σ, σ)
                    else
                        push!(y_A_Z_TKE.A, A)
                        push!(y_A_Z_TKE.Z, Z)
                        push!(y_A_Z_TKE.TKE, TKE)
                        push!(y_A_Z_TKE.Value, 2*val)
                        push!(y_A_Z_TKE.σ, σ*sqrt(2))
                    end
                end
            end
        end
    end
    for A in sort(unique(dY.A))
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
            for TKE in sort(unique(dY.TKE[dY.A .== A]))
                val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                if A != A_0 - A
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                elseif !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                end
            end
        end
    end
    #Renormalizes data
    f = 200/sum(y_A_Z_TKE.Value)
    y_A_Z_TKE.Value .= y_A_Z_TKE.Value .* f
    y_A_Z_TKE.σ .= y_A_Z_TKE.σ .* f
    return y_A_Z_TKE
end
#Average q(A,Z,TKE) over p(Z,A) so it becomes q(A,TKE)
function Average_over_Z(q_A_Z_TKE::Distribution, fragmdomain::Distribution)
    q_A_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(q_A_Z_TKE.A)
        for TKE in unique(q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A)])
            Denominator = 0.0
            Numerator = 0.0
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
                    Denominator += P_A_Z
                    Numerator += P_A_Z * value
                end
            end
            if Denominator > 0
                push!(q_A_TKE.A, A)
                push!(q_A_TKE.TKE, TKE)
                push!(q_A_TKE.Value, Numerator/Denominator)
            end
        end
    end
    return q_A_TKE
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(A)
function Average_over_TKE_Z(q_A_Z_TKE::Distribution, y_A_Z_TKE::Distribution)
    q_A = Distribution_unidym(Int[], Float64[], Float64[])
    for A in unique(q_A_Z_TKE.A)
        Denominator = 0.0
        Numerator = 0.0
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE * value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_A.Argument, A)
            push!(q_A.Value, Numerator/Denominator)
        end
    end
    return q_A
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) so it becomes q(TKE)
function Average_over_A_Z(q_A_Z_TKE::Distribution, y_A_Z_TKE::Distribution)
    q_TKE = Distribution_unidym(Float64[], Float64[], Float64[])
    for TKE in unique(q_A_Z_TKE.TKE)
        Denominator = 0.0
        Numerator = 0.0
        for A in unique(q_A_Z_TKE.A[(q_A_Z_TKE.TKE .== TKE)])
            for Z in q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE * value
                    end
                end
            end
        end
        if Denominator > 0
            push!(q_TKE.Argument, TKE)
            push!(q_TKE.Value, Numerator/Denominator)
        end
    end
    return q_TKE
end
#Average q(A,Z,TKE) over Y(A,Z,TKE) to get average value <q>
function Average_value(q_A_Z_TKE::Distribution, y_A_Z_TKE::Distribution, mass_number_range)
    Denominator = 0.0
    Numerator = 0.0
    for A in mass_number_range
        for Z in unique(q_A_Z_TKE.Z[(q_A_Z_TKE.A .== A)])
            for TKE in q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.Z .== Z)]
                value = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== A) .& (q_A_Z_TKE.TKE .== TKE) .& (q_A_Z_TKE.Z .== Z)][1]
                if !isnan(value)
                    if isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        Denominator += Y_A_Z_TKE
                        Numerator += Y_A_Z_TKE * value
                    end
                end
            end
        end
    end
    return Numerator/Denominator
end
#Get L-H pair value from given q(A) distribution and A_H 
function Pair_value(q_A::Distribution_unidym, A_0, A_H)
    val_L = q_A.Value[q_A.Argument .== A_0-A_H][1]
    val_H = q_A.Value[q_A.Argument .== A_H][1]
    return val_L + val_H
end
#Obtain Y(Z,Aₚ,TKE), Y(Z,Aₚ) & Y(Aₚ) distributions from Y(A,Z,TKE) & ν(A,Z, TKE)
function Yield_post_neutron(y_A_Z_TKE::Distribution, ν_A_Z_TKE::Distribution)
    y_Aₚ_Z_TKE = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    y_Aₚ_Z = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    y_Aₚ = Distribution_unidym(Int[], Float64[], Float64[])
    for A in unique(y_A_Z_TKE.A)
        for Z in unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A])
            for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z)]
                Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                σY_A_Z_TKE = y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)][1]
                if isassigned(ν_A_Z_TKE.Value[(ν_A_Z_TKE.A .== A) .& (ν_A_Z_TKE.Z .== Z) .& (ν_A_Z_TKE.TKE .== TKE)], 1)
                    ν = ν_A_Z_TKE.Value[(ν_A_Z_TKE.A .== A) .& (ν_A_Z_TKE.Z .== Z) .& (ν_A_Z_TKE.TKE .== TKE)][1]
                    Aₚ = A - ν
                else 
                    Aₚ = A 
                end
                if !isassigned(y_Aₚ_Z_TKE.Value[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)], 1)
                    push!(y_Aₚ_Z_TKE.A, Aₚ)
                    push!(y_Aₚ_Z_TKE.Z, Z)
                    push!(y_Aₚ_Z_TKE.TKE, TKE)
                    push!(y_Aₚ_Z_TKE.Value, Y_A_Z_TKE)
                    push!(y_Aₚ_Z_TKE.σ, σY_A_Z_TKE)
                else
                    y_Aₚ_Z_TKE.Value[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)] .+= Y_A_Z_TKE
                    y_Aₚ_Z_TKE.σ[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)] .+= sqrt(sum(y_Aₚ_Z_TKE.σ[(y_Aₚ_Z_TKE.A .== Aₚ) .& (y_Aₚ_Z_TKE.Z .== Z) .& (y_Aₚ_Z_TKE.TKE .== TKE)].^2) + σY_A_Z_TKE^2)
                end
            end
        end
    end
    for A in sort(unique(y_Aₚ_Z_TKE.A))
        for Z in unique(y_Aₚ_Z_TKE.Z[(y_Aₚ_Z_TKE.A .== A)])
            Y_Aₚ_Z = sum(y_Aₚ_Z_TKE.Value[(y_Aₚ_Z_TKE.A .== A) .& (y_Aₚ_Z_TKE.Z .== Z)])
            σY_Aₚ_Z = sqrt(sum(y_Aₚ_Z_TKE.σ[(y_Aₚ_Z_TKE.A .== A) .& (y_Aₚ_Z_TKE.Z .== Z)].^2))
            push!(y_Aₚ_Z.A, A)
            push!(y_Aₚ_Z.Z, Z)
            push!(y_Aₚ_Z.Value, Y_Aₚ_Z)
            push!(y_Aₚ_Z.σ, σY_Aₚ_Z)
        end
        Y_Aₚ = sum(y_Aₚ_Z_TKE.Value[y_Aₚ_Z_TKE.A .== A])
        σY_Aₚ = sqrt(sum(y_Aₚ_Z_TKE.σ[y_Aₚ_Z_TKE.A .== A].^2))
        push!(y_Aₚ.Argument, A)
        push!(y_Aₚ.Value, Y_Aₚ)
        push!(y_Aₚ.σ, σY_Aₚ)
    end
    return y_Aₚ_Z_TKE, y_Aₚ_Z, y_Aₚ
end
#####
println("*averaging data over $yield_distribution_filename experimental Yield distribution")

y_A_Z_TKE = Process_yield_data(A₀, fragmdomain, dY)
#=
ν_A_Z_TKE = Neutron_multiplicity_A_Z_TKE(DataFrame(
A = Raw_output_datafile.A,
Z = Raw_output_datafile.Z,
TKE = Raw_output_datafile.TKE,
No_Sequence = Raw_output_datafile.No_Sequence
))
avg_ν_L = Average_value(ν_A_Z_TKE, y_A_Z_TKE, LF_range)
avg_ν_H = Average_value(ν_A_Z_TKE, y_A_Z_TKE, HF_range)
avg_ν = (avg_ν_L + avg_ν_H)/2

T_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
    A = Raw_output_datafile.A,
    Z = Raw_output_datafile.Z,
    TKE = Raw_output_datafile.TKE,
    No_Sequence = Raw_output_datafile.No_Sequence,
    Value = Raw_output_datafile.Tₖ
))
avg_T_L = Average_value(T_A_Z_TKE, y_A_Z_TKE, LF_range)
avg_T_H = Average_value(T_A_Z_TKE, y_A_Z_TKE, HF_range)
avg_T = (avg_T_L + avg_T_H)/2

avg_ε_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
    A = Raw_output_datafile.A,
    Z = Raw_output_datafile.Z,
    TKE = Raw_output_datafile.TKE,
    No_Sequence = Raw_output_datafile.No_Sequence,
    Value = Raw_output_datafile.Avg_εₖ
))
avg_ε_L = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, LF_range)
avg_ε_H = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, HF_range)
avg_ε = (avg_ε_L + avg_ε_H)/2

open("output_data/$(fissionant_nucleus_identifier)_Average_quantities.OUT", "w") do file
    write(file, "<ν_L> = $avg_ν_L\n<ν_H> = $avg_ν_H\n<ν> = $avg_ν\n\n")
    write(file, "<T_L> = $avg_T_L\n<T_H> = $avg_T_H\n<T> = $avg_T\n\n")
    write(file, "<ε_L> = $avg_ε_L\n<ε_H> = $avg_ε_H\n<ε> = $avg_ε\n\n")
end

ν_A_TKE = Average_over_Z(ν_A_Z_TKE, fragmdomain)
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_nu_A_TKE.OUT", 
    DataFrame(A = ν_A_TKE.A, TKE = ν_A_TKE.TKE, ν = ν_A_TKE.Value), 
    writeheader=true, newline="\r\n", delim=' '
)

ν_A = Average_over_TKE_Z(ν_A_Z_TKE, y_A_Z_TKE)
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_nu_A.OUT", 
    DataFrame(A = ν_A.Argument, ν = ν_A.Value), 
    writeheader=true, newline="\r\n", delim=' '
)

ν_A_Pair = [Pair_value(ν_A, A₀, A_H) for A_H in ν_A.Argument[ν_A.Argument .>= A_H_min]]
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_nu_AH_Pair.OUT", 
    DataFrame(A = ν_A.Argument[ν_A.Argument .>= A_H_min], ν_Pair = ν_A_Pair), 
    writeheader=true, newline="\r\n", delim=' '
)

ν_TKE = Average_over_A_Z(ν_A_Z_TKE, y_A_Z_TKE)
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_nu_TKE.OUT", 
    DataFrame(TKE = ν_TKE.Argument, ν = ν_TKE.Value), 
    writeheader=true, newline="\r\n", delim=' '
)

y_Ap_Z_TKE, y_Ap_Z, y_Ap = Yield_post_neutron(y_A_Z_TKE, ν_A_Z_TKE)
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_Y_Ap.OUT", 
    DataFrame(Aₚ = y_Ap.Argument, Y = y_Ap.Value, σ = y_Ap.σ), 
    writeheader=true, newline="\r\n", delim=' '
)
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_Y_Ap_Z.OUT", 
    DataFrame(Aₚ = y_Ap_Z.A, Z = y_Ap_Z.Z, Y = y_Ap_Z.Value, σ = y_Ap_Z.σ), 
    writeheader=true, newline="\r\n", delim=' '
)
=#