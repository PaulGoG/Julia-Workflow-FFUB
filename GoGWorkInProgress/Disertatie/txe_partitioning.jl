#=
Different methods for partitioning the Total Excitation Energy (TXE) between the 
light and heavy fission fragments:

1. Modelling at scission using extra-deformation energies of the fragments (difference between 
deformation energy at scission and deformation energy at total acceleration of the fragments)

2. Partitioning ratios are provided via parametrization by segments

3. Constant temperature ratio R_T = T_L/T_H for fragments at total acceleration in Fermi Gas regime
=#
function TXE_partitioning(A_0, Z_0, A_H_min, A_H_max, Eₙ, fragmdomain, dΔE_def, tkerange, density_parameter_type, density_parameter_datafile, dm)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    if !isnan(Sₙ[1])
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                Z_L = Z_0 - Z_H
                if isassigned(dΔE_def.Value[(dΔE_def.A .== A_H) .& (dΔE_def.Z .== Z_H)], 1) && isassigned(dΔE_def.Value[(dΔE_def.A .== A_L) .& (dΔE_def.Z .== Z_L)], 1)
                    ΔE_def_H = dΔE_def.Value[(dΔE_def.A .== A_H) .& (dΔE_def.Z .== Z_H)][1]
                    ΔE_def_L = dΔE_def.Value[(dΔE_def.A .== A_L) .& (dΔE_def.Z .== Z_L)][1]
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    if !isnan(Q[1])
                        a_L = density_parameter(density_parameter_type, A_L, Z_L, density_parameter_datafile)
                        a_H = density_parameter(density_parameter_type, A_H, Z_H, density_parameter_datafile)
                        if !isnan(a_L) && !isnan(a_H)
                            r = a_L/a_H
                            for TKE in tkerange
                                TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, Sₙ[1], Sₙ[2], Eₙ)[1]
                                if TXE > 0
                                    E_scission = TXE - (ΔE_def_L + ΔE_def_H)
                                    if E_scission > 0
                                        E_excit_L = E_scission/(1 + r) + ΔE_def_L
                                        E_excit_H = r * E_scission/(1 + r) + ΔE_def_H
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
function TXE_partitioning(A_0, Z_0, A_H_min, A_H_max, Eₙ, fragmdomain, Points, tkerange, dm)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    if !isnan(Sₙ[1])
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            Ratio = Segments_TXE_partitioning(Points, A_H)
            if !isnan(Ratio)
                for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                    Z_L = Z_0 - Z_H
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    if !isnan(Q[1])
                        for TKE in tkerange
                            TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, Sₙ[1], Sₙ[2], Eₙ)[1]
                            if TXE > 0 
                                E_excit_H = TXE * Ratio
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
function TXE_partitioning(A_0, Z_0, A_H_min, A_H_max, Eₙ, fragmdomain, tkerange, density_parameter_type, density_parameter_datafile, dm)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    if !isnan(Sₙ[1])
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                Z_L = Z_0 - Z_H
                Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                a_L = density_parameter(density_parameter_type, A_L, Z_L, density_parameter_datafile)
                a_H = density_parameter(density_parameter_type, A_H, Z_H, density_parameter_datafile)
                if !isnan(Q[1]) && !isnan(a_L) && !isnan(a_H)
                    r = a_L/a_H
                    Ratio = 1/(1 + r*R_T²)
                    for TKE in tkerange
                        TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, Sₙ[1], Sₙ[2], Eₙ)[1]
                        if TXE > 0 
                            E_excit_H = TXE * Ratio
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
        Sort_TXE_partitioning(E_excit, fragmdomain)
        return E_excit
    else 
        error("Neutron separation energy for fissionant nucleus ($(A_0),$(Z_0)) could not be calculated!")
    end
end
function TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, density_parameter_type, density_parameter_datafile, dm)
    if txe_partitioning_type == "MSCZ"
        E_excitation = TXE_partitioning(A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, density_parameter_type, density_parameter_datafile, dm)
    elseif txe_partitioning_type == "PARAM"
        E_excitation = TXE_partitioning(A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, dm)
    elseif txe_partitioning_type == "RT"
        E_excitation = TXE_partitioning(A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, tkerange, density_parameter_type, density_parameter_datafile, dm)
    end
    return E_excitation
end
function Sort_TXE_partitioning(E_excit, fragmdomain)
    aux_A = copy(E_excit.A)
    aux_Z = copy(E_excit.Z)
    aux_TKE = zeros(length(E_excit.TKE))
    aux_Value = zeros(length(E_excit.Value))
    aux_index = 1
    for index_fragmdomain in eachindex(fragmdomain.A)
        for TKE in tkerange
            if isassigned(E_excit.Value[(E_excit.A .== fragmdomain.A[index_fragmdomain]) .& (E_excit.Z .== fragmdomain.Z[index_fragmdomain]) .& (E_excit.TKE .== TKE)], 1)
                aux_A[aux_index] = fragmdomain.A[index_fragmdomain]
                aux_Z[aux_index] = fragmdomain.Z[index_fragmdomain]
                aux_Value[aux_index] = E_excit.Value[(E_excit.A .== fragmdomain.A[index_fragmdomain]) .& (E_excit.Z .== fragmdomain.Z[index_fragmdomain]) .& (E_excit.TKE .== TKE)][1]
                aux_TKE[aux_index] = E_excit.TKE[(E_excit.A .== fragmdomain.A[index_fragmdomain]) .& (E_excit.Z .== fragmdomain.Z[index_fragmdomain]) .& (E_excit.TKE .== TKE)][1]
                aux_index += 1
            end
        end
    end
    for index in eachindex(aux_Value)
        E_excit.A[index] = aux_A[index]
        E_excit.Z[index] = aux_Z[index]
        E_excit.TKE[index] = aux_TKE[index]
        E_excit.Value[index] = aux_Value[index]
    end
end
function Segments_TXE_partitioning(P, x)
    for index in axes(P, 1)
        if P[index][1] > x
            x₀ = P[index-1][1]
            y₀ = P[index-1][2]
            x₁ = P[index][1]
            y₁ = P[index][2]
            slope = (y₁ - y₀)/(x₁ - x₀)
            y = y₀ + slope*(x - x₀)
            return y
        end
    end
    return NaN
end