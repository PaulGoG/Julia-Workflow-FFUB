#=
Different methods for partitioning the Total Excitation Energy (TXE) between the 
light and heavy fission fragments:

1. Modelling at scission using extra-deformation energies of the fragments (difference between 
deformation energy at scission and deformation energy at total acceleration of the fragments)

2. Partitioning ratios are provided via parametrization by segments

3. Temperature ratio R_T = T_L/T_H for fragments at total acceleration in Fermi Gas regime provided via parametrization by segments
=#
#####
#Modelling at scission
function TXE_partitioning(A_0, Z_0, A_H_range, fission_type, E_incident, fragmdomain::Distribution, dΔE_def::DataFrame, tkerange, density_parameter_type, density_parameter_data, dm::DataFrame)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    E_CN = Compound_nucleus_energy(fission_type, A_0, Z_0, E_incident, dm)
    if !isnan(E_CN[1])
        for A_H in A_H_range
            A_L = A_0 - A_H
            for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                Z_L = Z_0 - Z_H
                if isassigned(dΔE_def.Value[(dΔE_def.A .== A_H) .& (dΔE_def.Z .== Z_H)], 1) && isassigned(dΔE_def.Value[(dΔE_def.A .== A_L) .& (dΔE_def.Z .== Z_L)], 1)
                    ΔE_def_H = dΔE_def.Value[(dΔE_def.A .== A_H) .& (dΔE_def.Z .== Z_H)][1]
                    ΔE_def_L = dΔE_def.Value[(dΔE_def.A .== A_L) .& (dΔE_def.Z .== Z_L)][1]
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    if !isnan(Q[1])
                        a_L = density_parameter(density_parameter_type, A_L, Z_L, density_parameter_data)
                        a_H = density_parameter(density_parameter_type, A_H, Z_H, density_parameter_data)
                        if !isnan(a_L) && !isnan(a_H)
                            r = a_L/a_H
                            for TKE in tkerange
                                TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, E_CN[1], E_CN[2])[1]
                                if TXE > 0
                                    E_scission = TXE - (ΔE_def_L + ΔE_def_H)
                                    if E_scission > 0
                                        E_excit_L = r *E_scission/(1 + r) + ΔE_def_L
                                        E_excit_H = E_scission/(1 + r) + ΔE_def_H
                                        if A_L != A_H
                                            push!(E_excit.A, A_H)
                                            push!(E_excit.Z, Z_H)
                                            push!(E_excit.TKE, TKE)
                                            push!(E_excit.Value, E_excit_H)
                                            push!(E_excit.A, A_L)
                                            push!(E_excit.Z, Z_L)
                                            push!(E_excit.TKE, TKE)
                                            push!(E_excit.Value, E_excit_L)
                                        elseif Z_L != Z_H
                                            if !isassigned(E_excit.Value[(E_excit.A .== A_H) .& (E_excit.Z .== Z_H) .& (E_excit.TKE .== TKE)], 1)
                                                push!(E_excit.A, A_H)
                                                push!(E_excit.Z, Z_H)
                                                push!(E_excit.TKE, TKE)
                                                push!(E_excit.Value, E_excit_H)
                                                push!(E_excit.A, A_L)
                                                push!(E_excit.Z, Z_L)
                                                push!(E_excit.TKE, TKE)
                                                push!(E_excit.Value, E_excit_L)
                                            end
                                        else
                                            push!(E_excit.A, A_H)
                                            push!(E_excit.Z, Z_H)
                                            push!(E_excit.TKE, TKE)
                                            push!(E_excit.Value, TXE/2)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        Sort_TXE_partitioning(E_excit, fragmdomain)
        return E_excit
    else 
        error("Neutron separation energy for fissionant nucleus ($(A_0),$(Z_0)) could not be calculated!")
    end
end
#Direct partitioning ratio by segments
function TXE_partitioning(A_0, Z_0, A_H_range, fission_type, E_incident, fragmdomain::Distribution, Points::Vector{Tuple{Int64, Float64}}, tkerange, dm::DataFrame)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    E_CN = Compound_nucleus_energy(fission_type, A_0, Z_0, E_incident, dm)
    if !isnan(E_CN[1])
        for A_H in A_H_range
            A_L = A_0 - A_H
            Ratio = Segments_TXE_partitioning(Points, A_H)
            if Ratio <= 1
                for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                    Z_L = Z_0 - Z_H
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    if !isnan(Q[1])
                        for TKE in tkerange
                            TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, E_CN[1], E_CN[2])[1]
                            if TXE > 0 
                                E_excit_H = TXE *Ratio
                                E_excit_L = TXE - E_excit_H
                                if A_L != A_H
                                    push!(E_excit.A, A_H)
                                    push!(E_excit.Z, Z_H)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, E_excit_H)
                                    push!(E_excit.A, A_L)
                                    push!(E_excit.Z, Z_L)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, E_excit_L)
                                elseif Z_L != Z_H
                                    if !isassigned(E_excit.Value[(E_excit.A .== A_H) .& (E_excit.Z .== Z_H) .& (E_excit.TKE .== TKE)], 1)
                                        push!(E_excit.A, A_H)
                                        push!(E_excit.Z, Z_H)
                                        push!(E_excit.TKE, TKE)
                                        push!(E_excit.Value, E_excit_H)
                                        push!(E_excit.A, A_L)
                                        push!(E_excit.Z, Z_L)
                                        push!(E_excit.TKE, TKE)
                                        push!(E_excit.Value, E_excit_L)
                                    end
                                else
                                    push!(E_excit.A, A_H)
                                    push!(E_excit.Z, Z_H)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, TXE/2)
                                end
                            end
                        end
                    end
                end
            end
        end
        Sort_TXE_partitioning(E_excit, fragmdomain)
        return E_excit
    else 
        error("Neutron separation energy for fissionant nucleus ($(A_0),$(Z_0)) could not be calculated!")
    end
end
#R_T provided by segments
function TXE_partitioning(A_0, Z_0, A_H_range, fission_type, E_incident, fragmdomain::Distribution, Points::Vector{Tuple{Int64, Float64}}, tkerange, density_parameter_type, density_parameter_data, dm::DataFrame)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    E_CN = Compound_nucleus_energy(fission_type, A_0, Z_0, E_incident, dm)
    if !isnan(E_CN[1])
        for A_H in A_H_range
            A_L = A_0 - A_H
            R_T = Segments_TXE_partitioning(Points, A_H)
            if !isnan(R_T)
                for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                    Z_L = Z_0 - Z_H
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    a_L = density_parameter(density_parameter_type, A_L, Z_L, density_parameter_data)
                    a_H = density_parameter(density_parameter_type, A_H, Z_H, density_parameter_data)
                    if !isnan(Q[1]) && !isnan(a_L) && !isnan(a_H)
                        r = a_L/a_H
                        Ratio = 1/(1 + r *R_T^2)
                        for TKE in tkerange
                            TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, E_CN[1], E_CN[2])[1]
                            if TXE > 0 
                                E_excit_H = TXE *Ratio
                                E_excit_L = TXE - E_excit_H
                                if A_L != A_H
                                    push!(E_excit.A, A_H)
                                    push!(E_excit.Z, Z_H)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, E_excit_H)
                                    push!(E_excit.A, A_L)
                                    push!(E_excit.Z, Z_L)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, E_excit_L)
                                elseif Z_L != Z_H
                                    if !isassigned(E_excit.Value[(E_excit.A .== A_H) .& (E_excit.Z .== Z_H) .& (E_excit.TKE .== TKE)], 1)
                                        push!(E_excit.A, A_H)
                                        push!(E_excit.Z, Z_H)
                                        push!(E_excit.TKE, TKE)
                                        push!(E_excit.Value, E_excit_H)
                                        push!(E_excit.A, A_L)
                                        push!(E_excit.Z, Z_L)
                                        push!(E_excit.TKE, TKE)
                                        push!(E_excit.Value, E_excit_L)
                                    end
                                else
                                    push!(E_excit.A, A_H)
                                    push!(E_excit.Z, Z_H)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, TXE/2)
                                end
                            end
                        end
                    end
                end
            end
        end
        Sort_TXE_partitioning(E_excit, fragmdomain)
        return E_excit
    else 
        error("Neutron separation energy for fissionant nucleus ($(A_0),$(Z_0)) could not be calculated!")
    end
end
function TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_range, fission_type::String, E_incident, fragmdomain::Distribution, txe_partitioning_data, tkerange, density_parameter_type, density_parameter_data, dm::DataFrame)
    println("*partitioning Total Excitation Energy")
    if txe_partitioning_type == "MSCZ"
        return TXE_partitioning(A₀, Z₀, A_H_range, fission_type, E_incident, fragmdomain, txe_partitioning_data, tkerange, density_parameter_type, density_parameter_data, dm)
    elseif txe_partitioning_type == "PARAM"
        return TXE_partitioning(A₀, Z₀, A_H_range, fission_type, E_incident, fragmdomain, txe_partitioning_data, tkerange, dm)
    elseif txe_partitioning_type == "RT"
        return TXE_partitioning(A₀, Z₀, A_H_range, fission_type, E_incident, fragmdomain, txe_partitioning_data, tkerange, density_parameter_type, density_parameter_data, dm)
    end
end
function Sort_TXE_partitioning(E_excit, fragmdomain)
    aux_A = copy(E_excit.A)
    aux_Z = copy(E_excit.Z)
    aux_TKE = zeros(length(E_excit.TKE))
    aux_Value = zeros(length(E_excit.Value))
    aux_index = 1
    for index_fragmdomain in eachindex(fragmdomain.A)
        for TKE in unique(E_excit.TKE[(E_excit.A .== fragmdomain.A[index_fragmdomain]) .& (E_excit.Z .== fragmdomain.Z[index_fragmdomain])])
            aux_A[aux_index] = fragmdomain.A[index_fragmdomain]
            aux_Z[aux_index] = fragmdomain.Z[index_fragmdomain]
            aux_TKE[aux_index] = TKE
            aux_Value[aux_index] = E_excit.Value[(E_excit.A .== fragmdomain.A[index_fragmdomain]) .& (E_excit.Z .== fragmdomain.Z[index_fragmdomain]) .& (E_excit.TKE .== TKE)][1]
            aux_index += 1
        end
    end
    for index in eachindex(aux_Value)
        E_excit.A[index] = aux_A[index]
        E_excit.Z[index] = aux_Z[index]
        E_excit.TKE[index] = aux_TKE[index]
        E_excit.Value[index] = aux_Value[index]
    end
end
function Segments_TXE_partitioning(P::Vector{Tuple{Int, Float64}}, x::Int)
    for index in axes(P, 1)
        if P[index][1] >= x && index != 1
            x₀ = P[index-1][1]
            y₀ = P[index-1][2]
            x₁ = P[index][1]
            y₁ = P[index][2]
            slope = (y₁ - y₀)/(x₁ - x₀)
            y = y₀ + slope *(x - x₀)
            if y > 0
                return y
            else
                return NaN
            end
        end
    end
    return NaN
end