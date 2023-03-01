#=
Different methods for partitioning the Total Excitation Energy (TXE) between the 
light and heavy fission fragments:

1. Modelling at scission using extra-deformation energies of the fragments (difference between 
deformation energy at scission and deformation energy at total acceleration of the fragments)

2. Partitioning ratios are directly provided via datafile

3. Constant temperature ratio R_T = T_L/T_H for fragments at total acceleration in Fermi Gas regime
=#
function TXE_partitioning(A_0, Z_0, A_H_min, A_H_max, Eₙ, fragmdomain, dΔE_def, tkerange, density_parameter_type, density_parameter_datafile, dm)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    if !isnan(Sₙ[1])
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for index_Z_H in eachindex(fragmdomain.Z[fragmdomain.A .== A_H])
                Z_H = fragmdomain.Z[fragmdomain.A .== A_H][index_Z_H]
                Z_L = Z_0 - Z_H
                if isassigned(dΔE_def.Value[(dΔE_def.A .== A_H) .& (dΔE_def.Z .== Z_H)], 1) && isassigned(dΔE_def.Value[(dΔE_def.A .== A_L) .& (dΔE_def.Z .== Z_L)], 1)
                    ΔE_def_H = dΔE_def.Value[(dΔE_def.A .== A_H) .& (dΔE_def.Z .== Z_H)][1]
                    ΔE_def_L = dΔE_def.Value[(dΔE_def.A .== A_L) .& (dΔE_def.Z .== Z_L)][1]
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    Sₙ_L = Separation_energy(1, 0, A_L, Z_L, dm)[1]
                    Sₙ_H = Separation_energy(1, 0, A_H, Z_H, dm)[1]
                    if !isnan(Q[1]) && !isnan(Sₙ_L) && !isnan(Sₙ_H)
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
                                        if E_excit_H > Sₙ_H
                                            push!(E_excit.A, A_H)
                                            push!(E_excit.Z, Z_H)
                                            push!(E_excit.TKE, TKE)
                                            push!(E_excit.Value, E_excit_H)
                                        end
                                        if A_L != A_H_min && E_excit_L > Sₙ_L
                                            push!(E_excit.A, A_L)
                                            push!(E_excit.Z, Z_L)
                                            push!(E_excit.TKE, TKE)
                                            push!(E_excit.Value, E_excit_L)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return E_excit
end
function TXE_partitioning(A_0, Z_0, A_H_min, A_H_max, Eₙ, fragmdomain, dRatio, tkerange, dm)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    if !isnan(Sₙ[1])
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for index_Z_H in eachindex(fragmdomain.Z[fragmdomain.A .== A_H])
                Z_H = fragmdomain.Z[fragmdomain.A .== A_H][index_Z_H]
                Z_L = Z_0 - Z_H
                if isassigned(dRatio.Value[(dRatio.A .== A_H) .& (dRatio.Z .== Z_H)], 1)
                    Ratio = dRatio.Value[(dRatio.A .== A_H) .& (dRatio.Z .== Z_H)][1]
                    Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                    Sₙ_L = Separation_energy(1, 0, A_L, Z_L, dm)[1]
                    Sₙ_H = Separation_energy(1, 0, A_H, Z_H, dm)[1]
                    if !isnan(Q[1]) && !isnan(Sₙ_L) && !isnan(Sₙ_H)
                        for TKE in tkerange
                            TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, Sₙ[1], Sₙ[2], Eₙ)[1]
                            if TXE > 0 
                                E_excit_H = TXE * Ratio
                                E_excit_L = TXE - E_excit_H
                                if E_excit_H > Sₙ_H
                                    push!(E_excit.A, A_H)
                                    push!(E_excit.Z, Z_H)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, E_excit_H)
                                end
                                if A_L != A_H_min && E_excit_L > Sₙ_L
                                    push!(E_excit.A, A_L)
                                    push!(E_excit.Z, Z_L)
                                    push!(E_excit.TKE, TKE)
                                    push!(E_excit.Value, E_excit_L)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return E_excit
end
function TXE_partitioning(A_0, Z_0, A_H_min, A_H_max, Eₙ, fragmdomain, tkerange, density_parameter_type, density_parameter_datafile, dm)
    E_excit = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    if !isnan(Sₙ[1])
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for index_Z_H in eachindex(fragmdomain.Z[fragmdomain.A .== A_H])
                Z_H = fragmdomain.Z[fragmdomain.A .== A_H][index_Z_H]
                Z_L = Z_0 - Z_H
                Q = Q_value_released(A_0, Z_0, A_H, Z_H, dm)
                Sₙ_L = Separation_energy(1, 0, A_L, Z_L, dm)[1]
                Sₙ_H = Separation_energy(1, 0, A_H, Z_H, dm)[1]
                a_L = density_parameter(density_parameter_type, A_L, Z_L, density_parameter_datafile)
                a_H = density_parameter(density_parameter_type, A_H, Z_H, density_parameter_datafile)
                if !isnan(Q[1]) && !isnan(Sₙ_L) && !isnan(Sₙ_H) && !isnan(a_L) && !isnan(a_H)
                    r = a_L/a_H
                    Ratio = 1/(1 + r*R_T²)
                    for TKE in tkerange
                        TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, Sₙ[1], Sₙ[2], Eₙ)[1]
                        if TXE > 0 
                            E_excit_H = TXE * Ratio
                            E_excit_L = TXE - E_excit_H
                            if E_excit_H > Sₙ_H
                                push!(E_excit.A, A_H)
                                push!(E_excit.Z, Z_H)
                                push!(E_excit.TKE, TKE)
                                push!(E_excit.Value, E_excit_H)
                            end
                            if A_L != A_H_min && E_excit_L > Sₙ_L
                                push!(E_excit.A, A_L)
                                push!(E_excit.Z, Z_L)
                                push!(E_excit.TKE, TKE)
                                push!(E_excit.Value, E_excit_L)
                            end
                        end
                    end
                end
            end
        end
    end
    return E_excit
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