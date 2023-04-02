#=
Function bodies for solving the DSE energy conservation equations coresponding to 
constant and variable neutron evaporation cross section types.
=#
#####
#CONSTANT σₙ function
function DSE_equation_solver_CONSTANT_cs(E_excitation::Distribution, density_parameter_type, density_parameter_data, dm::DataFrame)
    Tₖ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    aₖ = Float64[]
    for A in unique(E_excitation.A)
        for Z in unique(E_excitation.Z[E_excitation.A .== A])
            Sₙ = Separation_energy(1, 0, A, Z, dm)[1]
            a_1 = density_parameter(density_parameter_type, A - 1, Z, density_parameter_data)
            if !isnan(a_1)
                for TKE in E_excitation.TKE[(E_excitation.A .== A) .& (E_excitation.Z .== Z)]
                    Eᵣ_k_last = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last = Sₙ
                    a_k = a_1
                    k = 1
                    Sum_avg_ε = 0.0
                    while Eᵣ_k_last - Sₙ_k_last > 0
                        T_k = (sqrt(1 + a_k *(Eᵣ_k_last - Sₙ_k_last)) - 1) /a_k
                        Sum_avg_ε += Average_neutron_energy(T_k)
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, T_k)
                        push!(aₖ, a_k)
                        push!(Tₖ.No_Sequence, k)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last = Energy_FermiGas(a_k, T_k)
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_data)
                    end
                    Eᵣ_k_last += Sum_avg_ε
                    while Eᵣ_k_last - Sₙ_k_last > 0
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, NaN)
                        push!(aₖ, a_k)
                        push!(Tₖ.No_Sequence, k)
                        Eᵣ_k_last -= Sₙ_k_last
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_data)
                    end
                    if k == 1
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, NaN)
                        push!(aₖ, a_k)
                        push!(Tₖ.No_Sequence, 0)
                    end
                end  
            end          
        end
    end
    return Tₖ, aₖ
end
#VARIABLE σₙ functions
#Parametrization for the force function S₀ of the s-wave neutron
function Neutron_strength_function(A)
    if A <= 23
        return 2.5e-5
    elseif A > 23 && A <= 55
        return 1.29695e-5 *A - 2.7328e-4
    elseif A > 55 && A <= 69
        return -2.1429e-5 *A + 1.6186e-3
    elseif A > 69 && A <= 79
        return 1.4e-4
    elseif A > 79 && A <= 90
        return -7.8182e-6 *A + 7.5764e-4
    elseif A > 90 && A <= 110
        return 5.4e-5
    elseif A > 110 && A <= 120
        return -4.4e-6 *A + 5.38e-4 
    elseif A > 120 && A <= 135
        return 6.6667e-6 *A - 7.9e-4
    elseif A > 135 && A <= 146
        return 2.4545e-5 *A - 3.2036e-3
    elseif A > 146 && A <= 152
        return 3.8e-4
    elseif A > 152 && A <= 158
        return -3.8333e-5 *A + 6.2067e-3
    elseif A > 158 && A <= 173
        return 1.5e-4
    elseif A > 173 && A <= 186
        return 1.0769e-5 *A - 1.7131e-3
    elseif A > 186 && A <= 208
        return -8.6364e-6 *A + 1.8964e-3
    elseif A > 208
        return 1e-4
    end
end
#Solves the transcendental equation for given sequence
function Solve_transcendental_eq(Eᵣ_k_last, Sₙ_k_last, a_k, A, k)
    σ₀ = π *(r₀ *0.1)^2 *(A - k)^(2/3)
    S₀ = Neutron_strength_function(A - k)
    αₖ = 10 *C_α *S₀ /σ₀
    epsilon_SCM(Tₖ) = Tₖ*(2*sqrt(Tₖ) + αₖ*3*sqrt(π)/4)/(sqrt(Tₖ) + αₖ*sqrt(π)/2)
    trans_eq(Tₖ) = a_k *Tₖ^2 + epsilon_SCM(Tₖ) + (Sₙ_k_last - Eᵣ_k_last)
    T_k = find_zero(trans_eq, (0.0, Inf))
    #TEST!
    if T_k > 2
        T_k = NaN
    end
    #TEST!
    return T_k, αₖ
end
function DSE_equation_solver_VARIABLE_cs(E_excitation::Distribution, density_parameter_type, density_parameter_data, dm::DataFrame)
    Tₖ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    aₖ = Float64[]
    αₖ = Float64[]
    for A in unique(E_excitation.A)
        for Z in unique(E_excitation.Z[E_excitation.A .== A])
            Sₙ = Separation_energy(1, 0, A, Z, dm)[1]
            a_1 = density_parameter(density_parameter_type, A - 1, Z, density_parameter_data)
            if !isnan(a_1)
                for TKE in E_excitation.TKE[(E_excitation.A .== A) .& (E_excitation.Z .== Z)]
                    Eᵣ_k_last = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last = Sₙ
                    a_k = a_1
                    k = 1
                    Sum_avg_ε = 0.0
                    while Eᵣ_k_last - Sₙ_k_last > 0
                        T_k, α_k = Solve_transcendental_eq(Eᵣ_k_last, Sₙ_k_last, a_k, A, k)
                        if !isnan(T_k)
                        Sum_avg_ε += Average_neutron_energy(α_k, T_k)
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, T_k)
                        push!(aₖ, a_k)
                        push!(αₖ, α_k)
                        push!(Tₖ.No_Sequence, k)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last = Energy_FermiGas(a_k, T_k)
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_data)
                        else
                            break
                        end
                    end
                    Eᵣ_k_last += Sum_avg_ε
                    while Eᵣ_k_last - Sₙ_k_last > 0
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, NaN)
                        push!(aₖ, a_k)
                        push!(αₖ, NaN)
                        push!(Tₖ.No_Sequence, k)
                        Eᵣ_k_last -= Sₙ_k_last
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_data)
                    end
                    if k == 1
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, NaN)
                        push!(aₖ, a_k)
                        push!(αₖ, NaN)
                        push!(Tₖ.No_Sequence, 0)
                    end
                end  
            end          
        end
    end
    return Tₖ, aₖ, αₖ
end
function DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_data, dm)
    println("*solving DSE energy conservation equations")
    if evaporation_cs_type == "CONSTANT"
        return DSE_equation_solver_CONSTANT_cs(E_excitation, density_parameter_type, density_parameter_data, dm)
    elseif evaporation_cs_type == "VARIABLE"
        return DSE_equation_solver_VARIABLE_cs(E_excitation, density_parameter_type, density_parameter_data, dm)
    end
end
