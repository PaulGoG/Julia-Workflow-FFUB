#=
Function bodies for solving the DSE conservation equations coresponding to 
constant and variable neutron evaporation cross section types.
=#
#####
function DSE_equation_solver_CONSTANT_cs(A_0, Z_0, A_H_min, A_H_max, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    Tₖ_L = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Tₖ_H = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[]) 
    aₖ_L = Float64[] 
    aₖ_H = Float64[]
    for A_H in A_H_min:A_H_max
        println("Progress: $(Int(round(100*(A_H - A_H_min + 1)/(A_H_max - A_H_min + 1))))%")
        A_L = A_0 - A_H
        for index_Z_H in eachindex(E_excitation.Z[E_excitation.A .== A_H])
            Z_H = E_excitation.Z[E_excitation.A .== A_H][index_Z_H]
            Z_L = Z_0 - Z_H
            Sₙ_L = Separation_energy(1, 0, A_L, Z_L, dm)[1]
            Sₙ_H = Separation_energy(1, 0, A_H, Z_H, dm)[1]
            a_1_H = density_parameter(density_parameter_type, A_0, Z_0, A_H - 1, Z_H, density_parameter_datafile)[2]
            a_1_L = density_parameter(density_parameter_type, A_0, Z_0, A_H + 1, Z_H, density_parameter_datafile)[1]
            for TKE in tkerange
                if isassigned(E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)], 1)
                    Eᵣ_k_last_H = E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last_H = Sₙ_H
                    a_k_H = a_1_H
                    k_H = 1
                    while Eᵣ_k_last_H > Sₙ_k_last_H && !isnan(a_k_H)
                        T_k_H = (sqrt(1 + a_k_H * (Eᵣ_k_last_H - Sₙ_k_last_H)) - 1)/a_k_H
                        push!(Tₖ_H.A, A_H)
                        push!(Tₖ_H.Z, Z_H)
                        push!(Tₖ_H.TKE, TKE)
                        push!(Tₖ_H.Value, T_k_H); 
                        push!(aₖ_H, a_k_H)
                        push!(Tₖ_H.NoSeq, k_H)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last_H = a_k_H * T_k_H^2
                        Sₙ_k_last_H = Separation_energy(1, 0, A_H - k_H, Z_H, dm)[1]
                        k_H += 1
                        a_k_H = density_parameter(density_parameter_type, A_0, Z_0, A_H - k_H, Z_H, density_parameter_datafile)[2]
                    end
                end
                if A_L != A_H
                    if isassigned(E_excitation.Value[(E_excitation.A .== A_L) .& (E_excitation.Z .== Z_L) .& (E_excitation.TKE .== TKE)], 1)
                        Eᵣ_k_last_L = E_excitation.Value[(E_excitation.A .== A_L) .& (E_excitation.Z .== Z_L) .& (E_excitation.TKE .== TKE)][1]
                        Sₙ_k_last_L = Sₙ_L
                        a_k_L = a_1_L
                        k_L = 1
                        while Eᵣ_k_last_L > Sₙ_k_last_L && !isnan(a_k_L)
                            T_k_L = (sqrt(1 + a_k_L * (Eᵣ_k_last_L - Sₙ_k_last_L)) - 1)/a_k_L
                            push!(Tₖ_L.A, A_L)
                            push!(Tₖ_L.Z, Z_L)
                            push!(Tₖ_L.TKE, TKE)
                            push!(Tₖ_L.Value, T_k_L);
                            push!(aₖ_L, a_k_L)
                            push!(Tₖ_L.NoSeq, k_L)
                            #Advance the sequence one step forward to be verified by the while loop
                            Eᵣ_k_last_L = a_k_L * T_k_L^2
                            Sₙ_k_last_L = Separation_energy(1, 0, A_L - k_L, Z_H, dm)[1]
                            k_L += 1
                            a_k_L = density_parameter(density_parameter_type, A_0, Z_0, A_H + k_L, Z_H, density_parameter_datafile)[2]
                        end
                    end
                end
            end            
        end
    end
    println("DSE energy conservation equations done!")
    return Tₖ_L, Tₖ_H, aₖ_L, aₖ_H
end
function DSE_equation_solver_VARIABLE_cs(A_0, Z_0, A_H_min, A_H_max, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    Tₖ_L = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    Tₖ_H = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[]) 
    aₖ_L = Float64[]
    aₖ_H = Float64[]
    αₖ_L = Float64[]
    αₖ_H = Float64[]
    for A_H in A_H_min:A_H_max
        println("Progress: $(Int(round(100*(A_H - A_H_min + 1)/(A_H_max - A_H_min + 1))))%")
        A_L = A_0 - A_H
        for index_Z_H in eachindex(E_excitation.Z[E_excitation.A .== A_H])
            Z_H = E_excitation.Z[E_excitation.A .== A_H][index_Z_H]
            Z_L = Z_0 - Z_H
            Sₙ_L = Separation_energy(1, 0, A_L, Z_L, dm)[1]
            Sₙ_H = Separation_energy(1, 0, A_H, Z_H, dm)[1]
            a_1_H = density_parameter(density_parameter_type, A_0, Z_0, A_H - 1, Z_H, density_parameter_datafile)[2]
            a_1_L = density_parameter(density_parameter_type, A_0, Z_0, A_H + 1, Z_H, density_parameter_datafile)[1]
            for TKE in tkerange
                if isassigned(E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)], 1)
                    Eᵣ_k_last_H = E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last_H = Sₙ_H
                    a_k_H = a_1_H
                    k_H = 1
                    while Eᵣ_k_last_H > Sₙ_k_last_H && !isnan(a_k_H)
                        T_k_H, α_k_H = Solve_transcendental_eq(Eᵣ_k_last_H, Sₙ_k_last_H, a_k_H, A_H, k_H)
                        push!(Tₖ_H.A, A_H)
                        push!(Tₖ_H.Z, Z_H)
                        push!(Tₖ_H.TKE, TKE)
                        push!(Tₖ_H.Value, T_k_H)
                        push!(aₖ_H, a_k_H)
                        push!(αₖ_H, α_k_H)
                        push!(Tₖ_H.NoSeq, k_H)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last_H = a_k_H * T_k_H^2
                        Sₙ_k_last_H = Separation_energy(1, 0, A_H - k_H, Z_H, dm)[1]
                        k_H += 1
                        a_k_H = density_parameter(density_parameter_type, A_0, Z_0, A_H - k_H, Z_H, density_parameter_datafile)[2]
                    end
                end
                if A_L != A_H
                    if isassigned(E_excitation.Value[(E_excitation.A .== A_L) .& (E_excitation.Z .== Z_L) .& (E_excitation.TKE .== TKE)], 1)
                        Eᵣ_k_last_L = E_excitation.Value[(E_excitation.A .== A_L) .& (E_excitation.Z .== Z_L) .& (E_excitation.TKE .== TKE)][1]
                        Sₙ_k_last_L = Sₙ_L
                        a_k_L = a_1_L
                        k_L = 1
                        while Eᵣ_k_last_L > Sₙ_k_last_L && !isnan(a_k_L)
                            T_k_L, α_k_L = Solve_transcendental_eq(Eᵣ_k_last_L, Sₙ_k_last_L, a_k_L, A_L, k_L)
                            push!(Tₖ_L.A, A_L)
                            push!(Tₖ_L.Z, Z_L)
                            push!(Tₖ_L.TKE, TKE)
                            push!(Tₖ_L.Value, T_k_L)
                            push!(aₖ_L, a_k_L)
                            push!(αₖ_L, α_k_L)
                            push!(Tₖ_L.NoSeq, k_L)
                            #Advance the sequence one step forward to be verified by the while loop
                            Eᵣ_k_last_L = a_k_L * T_k_L^2
                            Sₙ_k_last_L = Separation_energy(1, 0, A_L - k_L, Z_H, dm)[1]
                            k_L += 1
                            a_k_L = density_parameter(density_parameter_type, A_0, Z_0, A_H + k_L, Z_H, density_parameter_datafile)[2]
                        end
                    end
                end
            end            
        end
    end
    println("DSE energy conservation equations done!")
    return Tₖ_L, Tₖ_H, aₖ_L, aₖ_H
end
function DSE_equation_solver(evaporation_cs_type, A_0, Z_0, A_H_min, A_H_max, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    if evaporation_cs_type == "CONSTANT"
        DSE_Output = DSE_equation_solver_CONSTANT_cs(A_0, Z_0, A_H_min, A_H_max, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    elseif evaporation_cs_type == "VARIABLE"
        DSE_Output = DSE_equation_solver_VARIABLE_cs(A_0, Z_0, A_H_min, A_H_max, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    end
    return DSE_Output
end